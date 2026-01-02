function my_dipplotWithDensity(varargin)

% Default options
options = struct('dipxyz_clust_table', [], 'groupRanges', [], 'numGroups', [], ...
                 'clusters', [], 'color', [], 'plot', 'off', 'mri', '', ...
                 'view', 1, 'methodparam', 10, 'cmin_cmax', [], 'save', 0, ...
                 'cluster_data_all', [], 'rgb_matrix_all', []);

% Parse varargin as name-value pairs
for n = 1:2:length(varargin)
    if ischar(varargin{n})
        options.(varargin{n}) = varargin{n+1};
    end
end

% Extract options
dipxyz_clust_table = options.dipxyz_clust_table;
groupRanges = options.groupRanges;
numGroups = options.numGroups;
cluster_indices = options.clusters;
colors = options.color;
plot_opt = options.plot;
mri_file = options.mri;
view_opt = options.view;
methodparam = options.methodparam;
cmin_cmax = options.cmin_cmax;
save_opt = options.save;
cluster_data_all = options.cluster_data_all;
rgb_matrix_all = options.rgb_matrix_all;

% Validate inputs
if isempty(dipxyz_clust_table) || ~istable(dipxyz_clust_table)
    error('dipxyz_clust_table must be a non-empty table.');
end
if isempty(groupRanges) || length(groupRanges) < numGroups
    error('groupRanges must be a cell array with at least numGroups elements.');
end
if isempty(numGroups) || numGroups < 1
    error('numGroups must be a positive integer.');
end
if isempty(cluster_indices) || ~isnumeric(cluster_indices)
    error('Invalid or empty cluster indices.');
end

% Define the MRI file path
if isempty(mri_file)
    mri_file = 'C:\Users\docko\Desktop\UH\lab\eeg\eeg preprocessing\eeglab2024.1\plugins\dipfit5.5\standard_BEM\standard_mri.mat';
end

% Check if the MRI file exists and load it
if ~exist(mri_file, 'file')
    error('MRI file not found at: %s. Please check the path.', mri_file);
end
fprintf('MRI file found: %s\n', mri_file);
mri = load(mri_file); % Load the MRI structure
% Determine the correct field for MRI data
mri_fields = fieldnames(mri);
if ismember('mri', mri_fields) && isstruct(mri.mri)
    fprintf('Subfields in mri.mri: %s\n', strjoin(fieldnames(mri.mri), ', '));
    mri_data = mri.mri; % Use the full mri.mri structure for dipoledensity
    if isfield(mri.mri, 'anatomy')
        mri_volume = mri.mri.anatomy; % 3D MRI data for brainBlobBrowser
        fprintf('MRI anatomy dimensions: %s\n', mat2str(size(mri_volume)));
    else
        fprintf('Warning: No "anatomy" field in mri.mri. brainBlobBrowser may need adjustment.\n');
        mri_volume = []; % Fallback
    end
else
    fprintf('Warning: No "mri" field in MRI file. Using file path instead.\n');
    mri_data = mri_file; % Fallback to file path
    mri_volume = mri_file; % Fallback for brainBlobBrowser
end

% Debug: Check mri_data type
fprintf('Type of mri_data: %s\n', class(mri_data));

% Initialize accumulators for combined plots
allCentroidsPos = [];
allCentroidsMom = [];
allCentroidsRv = [];
allDipolesPos = [];
allDipolesMom = [];
allDipolesRv = [];
totalDipoles = 0;

% Preallocate newDipoles as a cell array to handle variable dipole counts
numClusters = length(cluster_indices);
newDipoles = cell(numClusters, 1);

% Initialize allDipColors and allDipNames as cell arrays
allDipColors = {};
allDipNames = {};

% Process each cluster
for n = 1:length(cluster_indices)
    tmpCluster = cluster_indices(n);
    fprintf('Processing Cluster %d\n', tmpCluster);

    % Use precomputed cluster data if available
    try
        fprintf('Extracting data for Cluster %d...\n', tmpCluster);
        if ~isempty(cluster_data_all) && ~isempty(cluster_data_all{n})
            cluster_data = cluster_data_all{n};
            rgbMatrix = rgb_matrix_all{n};
            fprintf('Using precomputed data for Cluster %d: %d rows\n', tmpCluster, size(cluster_data, 1));
        else
            idx = dipxyz_clust_table.ClustersID == tmpCluster;
            cluster_data = dipxyz_clust_table(idx, 5:7);
            if isempty(cluster_data)
                fprintf('Warning: No data for Cluster %d. Skipping.\n', tmpCluster);
                continue;
            end
            fprintf('Computed new data for Cluster %d: %d rows\n', tmpCluster, height(cluster_data));
            cluster_mean = mean(cluster_data{:,:}, 1);
            cluster_data = [cluster_data; array2table(cluster_mean, 'VariableNames', {'posx', 'posy', 'posz'})];
            cluster_data = table2array(cluster_data);
            rgbMatrix = zeros(size(cluster_data, 1), 3);
            rgbMatrix(end, :) = [0 0 1]; % Blue centroid
        end
    catch e
        fprintf('Error extracting data for Cluster %d: %s\n', tmpCluster, e.message);
        disp(e);
        rethrow(e);
    end

    % Compute centroid and statistics
    try
        fprintf('Computing statistics for Cluster %d...\n', tmpCluster);
        coordinates = cluster_data(1:end-1, :); % Exclude centroid
        if ~isempty(coordinates)
            centroid = mean(coordinates, 1);
            standardDeviation = std(coordinates);
            standardError = standardDeviation / sqrt(size(coordinates, 1));
        else
            centroid = [0 0 0];
            standardDeviation = [0 0 0];
            standardError = [0 0 0];
        end
        clusterReport = sprintf('\nCluster: %d\nCentroid in MNI:    [%2.0f %2.0f %2.0f]\nStandard Deviation: [%2.0f %2.0f %2.0f]\nStandard Error:    [%2.0f %2.0f %2.0f]\n', ...
            tmpCluster, centroid(1), centroid(2), centroid(3), ...
            standardDeviation(1), standardDeviation(2), standardDeviation(3), ...
            standardError(1), standardError(2), standardError(3));
        fprintf('%s', clusterReport);
    catch e
        fprintf('Error computing statistics for Cluster %d: %s\n', tmpCluster, e.message);
        disp(e);
        rethrow(e);
    end

    % Accumulate centroids and dipoles for combined plots
    try
        fprintf('Accumulating centroids for Cluster %d...\n', tmpCluster);
        allCentroidsPos = [allCentroidsPos; centroid];
        allCentroidsMom = [allCentroidsMom; zeros(1, 3)];
        allCentroidsRv = [allCentroidsRv; 0];
        if ~isempty(coordinates)
            numDipoles = size(coordinates, 1);
            startIdx = totalDipoles + 1;
            endIdx = totalDipoles + numDipoles;
            allDipolesPos(startIdx:endIdx, :) = coordinates;
            allDipolesMom(startIdx:endIdx, :) = zeros(numDipoles, 3);
            allDipolesRv(startIdx:endIdx, :) = zeros(numDipoles, 1);
            totalDipoles = endIdx;
            allDipNames = [allDipNames arrayfun(@(x) sprintf('Cluster %d, Dipole %d', tmpCluster, x), 1:numDipoles, 'UniformOutput', false)];
            cluster_colors = lines(length(cluster_indices));
            cluster_color = cluster_colors(n, :);
            allDipColors = [allDipColors repmat({cluster_color}, 1, numDipoles)];
        end
    catch e
        fprintf('Error accumulating centroids/dipoles for Cluster %d: %s\n', tmpCluster, e.message);
        disp(e);
        rethrow(e);
    end

    % Validate coordinates
    try
        fprintf('Validating coordinates for Cluster %d...\n', tmpCluster);
        if isempty(coordinates) || ~isnumeric(coordinates) || any(isnan(coordinates(:)))
            fprintf('Warning: Coordinates for Cluster %d are invalid or empty. Skipping plot.\n', tmpCluster);
            continue;
        end
        fprintf('Coordinates for Cluster %d are valid: %d rows\n', tmpCluster, size(coordinates, 1));
    catch e
        fprintf('Error validating coordinates for Cluster %d: %s\n', tmpCluster, e.message);
        disp(e);
        rethrow(e);
    end

    % Plot dipole density if plotting is enabled
    if strcmp(plot_opt, 'on')
        % Use a simple colormap
        customJet = jet(128);
        plotargs = {'cmap', customJet};

        % Prepare dipole data structure
        try
            fprintf('Creating dipModel for Cluster %d...\n', tmpCluster);
            dipModel = struct('posxyz', num2cell(coordinates, 2), ...
                              'momxyz', {zeros(size(coordinates, 1), 3)}, ...
                              'rv', {zeros(size(coordinates, 1), 1)});
        catch e
            fprintf('Error creating dipModel for Cluster %d: %s\n', tmpCluster, e.message);
            disp(e);
            rethrow(e);
        end

        % Try plotting with detailed error catching
        try
            fprintf('Running dipoledensity for Cluster %d...\n', tmpCluster);
            [dens3d, ~] = dipoledensity(dipModel, 'coordformat', 'MNI', ...
                'methodparam', methodparam, 'plot', 'on', 'norm2JointProb', 'on', ...
                'plotargs', plotargs, 'mri', mri_data);
            fprintf('dipoledensity completed successfully for Cluster %d\n', tmpCluster);
        catch e
            fprintf('Error in dipoledensity for Cluster %d: %s\n', tmpCluster, e.message);
            disp(e);
            rethrow(e);
        end

        % Validate and preprocess density data for brainBlobBrowser
        try
            fprintf('Checking density data for Cluster %d...\n', tmpCluster);
            if isempty(dens3d) || ~iscell(dens3d) || isempty(dens3d{1})
                fprintf('Warning: No valid density data for Cluster %d. Skipping further plotting.\n', tmpCluster);
                continue;
            end
            densityData = dens3d{1};
            fprintf('Original dimensions of dens3d{1} for Cluster %d: %s\n', tmpCluster, mat2str(size(densityData)));

            % Resample densityData to [91 109 91] to match brainBlobBrowser expectation
            target_dims = [91 109 91];
            [X, Y, Z] = meshgrid(1:size(densityData, 2), 1:size(densityData, 1), 1:size(densityData, 3));
            [Xq, Yq, Zq] = meshgrid(linspace(1, size(densityData, 2), target_dims(2)), ...
                                  linspace(1, size(densityData, 1), target_dims(1)), ...
                                  linspace(1, size(densityData, 3), target_dims(3)));
            densityData_resampled = interp3(X, Y, Z, densityData, Xq, Yq, Zq, 'linear', 0);
            fprintf('Resampled dimensions of densityData for Cluster %d: %s\n', tmpCluster, mat2str(size(densityData_resampled)));

            % Normalize densityData to [0, 1] as per brainBlobBrowser
            fprintf('Min density before normalization: %.7f, Max density: %.7f\n', min(densityData_resampled(:)), max(densityData_resampled(:)));
            if max(densityData_resampled(:)) > 0
                densityData_resampled = densityData_resampled / max(densityData_resampled(:));
            else
                fprintf('Warning: Max density is 0 for Cluster %d. Setting densityData to zeros.\n', tmpCluster);
                densityData_resampled = zeros(size(densityData_resampled));
            end
            fprintf('Max dipole density for Cluster %d after normalization: %.7f\n', tmpCluster, max(densityData_resampled(:)));

            % Ensure mri_volume is compatible
            if isempty(mri_volume) || ~isequal(size(mri_volume), target_dims)
                fprintf('Warning: mri_volume dimensions (%s) do not match target (%s). Using resampled mri_data.\n', ...
                    mat2str(size(mri_volume)), mat2str(target_dims));
                mri_volume_resampled = interp3(X, Y, Z, double(mri_volume), Xq, Yq, Zq, 'linear', 0);
            else
                mri_volume_resampled = mri_volume;
            end

        catch e
            fprintf('Error processing density data for Cluster %d: %s\n', tmpCluster, e.message);
            disp(e);
            rethrow(e);
        end

        % Set figure title
        h1 = gcf;
        set(h1, 'Name', sprintf('Cls %d Dipole Density Plot', tmpCluster));

        % Try brainBlobBrowser with resampled data
        try
            fprintf('Running brainBlobBrowser for Cluster %d...\n', tmpCluster);
            brainBlobBrowser('data', densityData_resampled, 'clusterReport', clusterReport, ...
                           'clusterStd', standardDeviation, 'mri', mri_volume_resampled, 'color', 'red');
            fprintf('brainBlobBrowser completed for Cluster %d\n', tmpCluster);
        catch e
            fprintf('Error in brainBlobBrowser for Cluster %d: %s\n', tmpCluster, e.message);
            fprintf('Stack trace for brainBlobBrowser error:\n');
            disp(e.stack);
            fprintf('Skipping brainBlobBrowser for Cluster %d and continuing...\n', tmpCluster);
            continue;
        end
        h2 = gcf;
        set(h2, 'Name', sprintf('Cls %d BrainBlobBrowser', tmpCluster));

        % Save figures if requested
        if save_opt == 1
            print(h1, '-deps', sprintf('dipDensity_Cls%d', tmpCluster), '-loose');
            print(h2, '-deps', sprintf('brainBlob_Cls%d', tmpCluster), '-loose');
        end
    end
end

% Debug print to verify centroids and dipoles
fprintf('Number of centroids accumulated: %d\n', size(allCentroidsPos, 1));
fprintf('Number of dipoles accumulated: %d\n', size(allDipolesPos, 1));

% Build allCentroids structure for dipplot
allCentroids = struct('posxyz', {}, 'momxyz', {}, 'rv', {});
for i = 1:size(allCentroidsPos, 1)
    allCentroids(i).posxyz = allCentroidsPos(i, :);
    allCentroids(i).momxyz = allCentroidsMom(i, :);
    allCentroids(i).rv = allCentroidsRv(i);
end

% Build allDipoles structure for dipplot
allDipoles = struct('posxyz', {}, 'momxyz', {}, 'rv', {});
for i = 1:size(allDipolesPos, 1)
    allDipoles(i).posxyz = allDipolesPos(i, :);
    allDipoles(i).momxyz = allDipolesMom(i, :);
    allDipoles(i).rv = allDipolesRv(i);
end

% Debug structure contents
fprintf('Inspecting allCentroids structure:\n');
for i = 1:length(allCentroids)
    fprintf('Centroid %d: posxyz = %s, momxyz = %s, rv = %.4f\n', ...
        i, mat2str(allCentroids(i).posxyz), mat2str(allCentroids(i).momxyz), allCentroids(i).rv);
end
fprintf('Inspecting allDipoles structure (first 5 entries):\n');
for i = 1:min(5, length(allDipoles))
    fprintf('Dipole %d: posxyz = %s, momxyz = %s, rv = %.4f\n', ...
        i, mat2str(allDipoles(i).posxyz), mat2str(allDipoles(i).momxyz), allDipoles(i).rv);
end

% Plot centroids for all clusters
if strcmp(plot_opt, 'on') && ~isempty(allCentroids)
    try
        fprintf('Plotting centroids for all clusters...\n');
        figure;
        dipplot(allCentroids, 'mri', mri_file, 'coordformat', 'MNI', ...
                'dipolelength', 0, 'spheres', 'on', 'projlines', 'on', ...
                'view', [45 35], 'projimg', 'off', 'color', lines(length(cluster_indices)));
        h3 = gcf;
        set(h3, 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off', 'menu', 'none', 'NumberTitle', 'off', ...
            'Name', 'Cluster Centroids');
        if save_opt == 1
            print(h3, '-deps', 'centroids', '-loose');
        end
        fprintf('Centroid plot generated successfully.\n');
    catch e
        fprintf('Error plotting centroids: %s\n', e.message);
        disp(e);
    end
end

% Plot all dipoles
if strcmp(plot_opt, 'on') && ~isempty(allDipoles)
    try
        fprintf('Plotting all dipoles...\n');
        figure;
        dipplot(allDipoles, 'mri', mri_file, 'coordformat', 'MNI', ...
                'dipolelength', 0, 'spheres', 'on', 'projlines', 'off', ...
                'view', [45 35], 'projimg', 'off', 'color', [allDipColors{:}]);
        h4 = gcf;
        set(h4, 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off', 'menu', 'none', 'NumberTitle', 'off', ...
            'Name', 'All Cluster Dipoles');
        if save_opt == 1
            print(h4, '-deps', 'all_dipoles', '-loose');
        end
        fprintf('All dipoles plot generated successfully.\n');
    catch e
        fprintf('Error plotting all dipoles: %s\n', e.message);
        disp(e);
    end
end

end % End of function

% Helper function: computecentroid (unchanged)
function dipole = computecentroid(alldipoles)
    max_r = 0;
    len = length(alldipoles);
    dipole.posxyz = [0 0 0];
    dipole.momxyz = [0 0 0];
    dipole.rv = 0;
    count = 0;
    numNaN = 0;
    warningon = 1;
    for k = 1:len 
        if size(alldipoles(k).posxyz, 1) == 2
            if all(alldipoles(k).posxyz(2, :) == [0 0 0])
                alldipoles(k).posxyz(2, :) = [];
                alldipoles(k).momxyz(2, :) = [];
            end
        end
        if ~isempty(alldipoles(k).posxyz)
            dipole.posxyz = dipole.posxyz + mean(alldipoles(k).posxyz, 1);
            dipole.momxyz = dipole.momxyz + mean(alldipoles(k).momxyz, 1);
            if ~isnan(alldipoles(k).rv)
                dipole.rv = dipole.rv + alldipoles(k).rv;
            else
                numNaN = numNaN + 1;
            end
            count = count + 1;
        elseif warningon
            disp('Some components do not have dipole information');
            warningon = 0;
        end
    end
    if count > 0
        dipole.posxyz = dipole.posxyz / count;
        dipole.momxyz = dipole.momxyz / count;
        dipole.rv = dipole.rv / (count - numNaN);
    else
        dipole.posxyz = [0 0 0];
        dipole.momxyz = [0 0 0];
        dipole.rv = 0;
    end
    if isfield(alldipoles, 'maxr')
        dipole.maxr = alldipoles(1).max_r;
    end
end