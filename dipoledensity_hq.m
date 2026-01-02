function [prob3d, mri] = dipoledensity_hq(dipoles, varargin)
% DIPOLEDENSITY_HQ - Standalone high-quality dipole density plotting
%
% This is a standalone version that does NOT depend on dipplot.
% Uses standard_mri.mat directly and handles coordinate transformations internally.
%
% Usage:
%   >> [dens3d, mri] = dipoledensity_hq(dipoles, 'key', val, ...);
%
% Inputs:
%   dipoles - Can be either:
%             1. A structure array with fields: posxyz, momxyz, rv
%             2. An Nx3 matrix of [x, y, z] coordinates in MNI space
%             3. A 3xN matrix of [x; y; z] coordinates in MNI space
%
% Optional key-value pairs:
%   'subjind'        - Subject index for each dipole (required for entropy methods)
%   'method'         - 'alldistance' (default), 'distance', 'entropy', 'relentropy'
%   'methodparam'    - Smoothing kernel std (mm) or number of dipoles (default: 20)
%   'weight'         - Weight for each dipole (default: ones)
%   'smooth'         - Additional 3D smoothing (default: 0)
%   'nsessions'      - Number of sessions for normalization (default: 1)
%   'subsample'      - Subsampling factor (default: 2)
%   'plot'           - 'on' or 'off' (default: 'on' if no output, else 'off')
%   'normalization'  - 'on' or 'off' (default: 'on')
%   'norm2JointProb' - 'on' or 'off' (default: 'off')
%   'volmesh_fname'  - Filename for volume mesh cache (default: 'volmesh_hq.mat')
%   'plotargs'       - Cell array of arguments for mri3dplot
%
% Outputs:
%   prob3d - 3D density array
%   mri    - MRI structure used for plotting
%
% Example:
%   % From structure array
%   [dens3d, mri] = dipoledensity_hq(dipModelAll);
%   
%   % From coordinate matrix
%   coords = [x_vals, y_vals, z_vals];  % Nx3 matrix in MNI coordinates
%   [dens3d, mri] = dipoledensity_hq(coords);
%
% Note: Input coordinates MUST be in MNI space (millimeters)
%
% Created by: Komal K. Kukkar
% Date: December 10, 2025

prob3d = []; mri = [];

if nargin < 1
   help dipoledensity_hq
   return
end

% Parse input parameters
g = finputcheck(varargin, { 
    'subjind'        'integer'  []               [];
    'method'         'string'   { 'relentropy','entropy','distance','alldistance' } 'alldistance';
    'methodparam'    'real'     []               20; 
    'weight'         { 'real','cell' }  []       [];
    'smooth'         'real'     []               0;
    'nsessions'      'integer'  []               1;
    'subsample'      'integer'  []               2;
    'plotargs'       'cell'     []               {};
    'plot'           'string'   { 'on','off' }   fastif(nargout == 0, 'on', 'off');
    'normalization'  'string'   { 'on','off' }   'on';
    'volmesh_fname'  'string'   []               'volmesh_hq.mat';
    'norm2JointProb' 'string'   { 'on','off' }   'off'});

if ischar(g), error(g); end

if ~strcmpi(g.method, 'alldistance') && isempty(g.subjind)
    error('Subject indices are required for entropy and distance methods');
end
if ~iscell(g.weight), g.weight = { g.weight }; end

% ========================================================================
% STEP 1: Convert input dipoles to standard format (get MNI coordinates)
% ========================================================================
fprintf('Processing dipole input...\n');

if isstruct(dipoles)
    % Input is structure array with posxyz field
    numDipoles = length(dipoles);
    allx = zeros(1, numDipoles);
    ally = zeros(1, numDipoles);
    allz = zeros(1, numDipoles);
    
    for i = 1:numDipoles
        if isfield(dipoles(i), 'posxyz') && length(dipoles(i).posxyz) >= 3
            allx(i) = dipoles(i).posxyz(1);
            ally(i) = dipoles(i).posxyz(2);
            allz(i) = dipoles(i).posxyz(3);
        else
            error('Dipole structure must have posxyz field with [x,y,z] coordinates');
        end
    end
    
    fprintf('Extracted %d dipoles from structure array\n', numDipoles);
    
elseif isnumeric(dipoles)
    % Input is numeric matrix
    if size(dipoles, 1) == 3 && size(dipoles, 2) ~= 3
        % 3xN format - transpose
        dipoles = dipoles';
    end
    
    if size(dipoles, 2) ~= 3
        error('Dipole matrix must be Nx3 or 3xN with [x,y,z] coordinates');
    end
    
    allx = dipoles(:, 1)';
    ally = dipoles(:, 2)';
    allz = dipoles(:, 3)';
    numDipoles = size(dipoles, 1);
    
    fprintf('Extracted %d dipoles from coordinate matrix\n', numDipoles);
else
    error('Dipoles must be either a structure array or numeric matrix');
end

% Verify coordinates are in reasonable MNI range
if max(abs([allx, ally, allz])) > 200
    warning('Dipole coordinates seem very large (>200mm). Are you sure they are in MNI space?');
end

fprintf('Dipole coordinate ranges:\n');
fprintf('  X: [%.1f, %.1f] mm\n', min(allx), max(allx));
fprintf('  Y: [%.1f, %.1f] mm\n', min(ally), max(ally));
fprintf('  Z: [%.1f, %.1f] mm\n', min(allz), max(allz));

% ========================================================================
% STEP 2: Set up weights and subject indices
% ========================================================================
if isempty(g.weight{1})
    g.weight = { ones(1, numDipoles) };
end

if ~iscell(g.weight)
    if length(g.weight) ~= numDipoles
        error('Weight vector must have %d elements (one per dipole)', numDipoles);
    end
else
    if length(g.weight{1}) ~= numDipoles || length(g.weight{end}) ~= numDipoles
        error('Weight vectors must have %d elements (one per dipole)', numDipoles);
    end
end

if isempty(g.subjind)
    g.subjind = ones(1, numDipoles);
else
    if length(g.subjind) ~= numDipoles
        error('Subject index vector must have %d elements (one per dipole)', numDipoles);
    end
end

% ========================================================================
% STEP 3: Load high-quality MRI (standard_mri.mat)
% ========================================================================
fprintf('\nLoading high-quality standard_mri.mat...\n');

% Find EEGLAB root
eeglabRoot = fileparts(which('eeglab'));
if isempty(eeglabRoot)
    error('EEGLAB not found in path. Please start EEGLAB first.');
end

% Look for standard_mri.mat in dipfit plugin
mriFile = fullfile(eeglabRoot, 'plugins', 'dipfit', 'standard_BEM', 'standard_mri.mat');

if ~exist(mriFile, 'file')
    fprintf('Standard MRI not found at expected location. Searching...\n');
    mriFiles = dir(fullfile(eeglabRoot, 'plugins', '**', 'standard_mri.mat'));
    if ~isempty(mriFiles)
        mriFile = fullfile(mriFiles(1).folder, mriFiles(1).name);
    else
        error('standard_mri.mat not found in EEGLAB dipfit plugin. Please ensure dipfit is installed.');
    end
end

fprintf('Using MRI: %s\n', mriFile);

% Load the MRI file
try
    mriLoaded = load('-mat', mriFile);
    if isfield(mriLoaded, 'mri')
        mri = mriLoaded.mri;
    else
        error('MRI structure not found in %s', mriFile);
    end
catch ME
    error('Failed to load standard_mri.mat: %s', ME.message);
end

fprintf('MRI loaded successfully\n');
fprintf('  Dimensions: [%d %d %d]\n', mri.dim);
fprintf('  Transform matrix:\n');
disp(mri.transform);

% ========================================================================
% STEP 4: Verify dipole coordinates are in correct space
% ========================================================================
% The dipoles should already be in MNI space (mm)
% standard_mri.mat uses this transform:
% [1  0  0  -91 ]
% [0  1  0  -126]
% [0  0  1  -73 ]
% [0  0  0   1  ]
%
% This means: MNI_coord = voxel_index * 1mm + offset
% So MNI coordinates should typically range from about -90 to +90 mm

fprintf('\nDipole coordinates are assumed to be in MNI space (millimeters)\n');
fprintf('No coordinate transformation will be applied.\n');

% ========================================================================
% STEP 5: Initialize density arrays
% ========================================================================
prob3d = {zeros(ceil(mri.dim/g.subsample)) };
for i = 2:length(g.weight)
    prob3d{i} = prob3d{1}; 
end

% Compute voxel volume
point1 = mri.transform * [ 1 1 1 1 ]';
point2 = mri.transform * [ 2 2 2 1 ]';
voxvol = sum((point1(1:3)-point2(1:3)).^2) * g.subsample^3; % in mm^3

fprintf('\nVoxel volume: %.2f mm^3\n', voxvol);

% ========================================================================
% STEP 6: Compute entropy parameters if needed
% ========================================================================
vals = unique_bc(g.subjind);
if strcmpi(g.method, 'relentropy') || strcmpi(g.method, 'entropy')
    newind = zeros(size(g.subjind));
    for index = 1:length(vals)
        tmpind = find(g.subjind == vals(index));
        totcount(index) = length(tmpind);
        newind(tmpind) = index;
    end
    g.subjind = newind;
    gp = totcount/sum(totcount);
    globent = -sum(gp.*log(gp));
    fprintf('Global entropy: %.4f\n', globent);
end

% ========================================================================
% STEP 7: Compute volume inside head mesh
% ========================================================================
fprintf('\nComputing brain volume mask...\n');

% Load standard BEM volume
dipfitdefs;
tmp = load('-mat', DIPOLEDENSITY_STDBEM);

if ~exist(g.volmesh_fname, 'file')
    if exist('ft_electroderealign', 'file') ~= 2
        error('dipoledensity_hq: FieldTrip toolbox is required'); 
    end
    
    fprintf('Computing volume mesh (this may take several minutes)...\n');
    
    % Create grid of voxel centers
    [X, Y, Z] = meshgrid(mri.xgrid(1:g.subsample:end) + g.subsample/2, ...
                         mri.ygrid(1:g.subsample:end) + g.subsample/2, ...
                         mri.zgrid(1:g.subsample:end) + g.subsample/2);
    [indX, indY, indZ] = meshgrid(1:length(mri.xgrid(1:g.subsample:end)), ...
                                  1:length(mri.ygrid(1:g.subsample:end)), ...
                                  1:length(mri.zgrid(1:g.subsample:end)));
    
    % Stack coordinates
    allpoints = [ X(:)'; Y(:)'; Z(:)' ];
    allinds   = [ indX(:)'; indY(:)'; indZ(:)' ];
    
    % Transform voxel indices to MNI coordinates
    allpoints = mri.transform * [ allpoints; ones(1, size(allpoints, 2)) ];
    allpoints(4,:) = [];
    
    fprintf('  Total voxels: %d\n', size(allpoints, 2));
    fprintf('  Coordinate range: X[%.1f,%.1f], Y[%.1f,%.1f], Z[%.1f,%.1f]\n', ...
        min(allpoints(1,:)), max(allpoints(1,:)), ...
        min(allpoints(2,:)), max(allpoints(2,:)), ...
        min(allpoints(3,:)), max(allpoints(3,:)));
    
    % Determine which voxels are inside the head
    olddir = pwd;
    tmppath = which('ft_electroderealign');
    tmppath = fullfile(fileparts(tmppath), 'private');
    cd(tmppath);
    inside = ft_inside_headmodel(allpoints', tmp.vol);
    Inside = find(inside); 
    Outside = find(~inside);
    cd(olddir);
    
    fprintf('  Voxels inside brain: %d (%.1f%%)\n', length(Inside), 100*length(Inside)/size(allpoints,2));
    
    % Save for future use
    try 
        save('-mat', g.volmesh_fname, 'allpoints', 'allinds', 'Inside', 'Outside');
        fprintf('  Saved volume mesh to: %s\n', g.volmesh_fname);
    catch
        warning('Could not save volume mesh file.');
    end
else
    fprintf('Loading pre-computed volume mesh from: %s\n', g.volmesh_fname);
    load('-mat', g.volmesh_fname);
    fprintf('  Voxels inside brain: %d\n', length(Inside));
end

InsidePoints  = allpoints(:, Inside);
InsideIndices = allinds(:, Inside);

% ========================================================================
% STEP 8: Compute density
% ========================================================================
fprintf('\nComputing dipole density...\n');
fprintf('Method: %s\n', g.method);
fprintf('Method parameter: %.2f\n', g.methodparam);

edges = [0.5:1:length(vals)+0.5];

if ~strcmpi(g.method, 'alldistance') 
    % Entropy or distance methods - scan voxels
    fprintf('Processing voxels (of %d): ', size(InsideIndices, 2));
    
    for i = 1:size(InsideIndices, 2)
        % Compute distances from this voxel to all dipoles
        alldists = (InsidePoints(1,i) - allx).^2 + ...
                   (InsidePoints(2,i) - ally).^2 + ...
                   (InsidePoints(3,i) - allz).^2;
        [tmpsort, indsort] = sort(alldists);
        tmpweights{1}   = g.weight{1}(indsort);
        tmpweights{end} = g.weight{end}(indsort);
       
        if strcmpi(g.method, 'relentropy') || strcmpi(g.method, 'entropy')
            % Entropy calculation
            subjs = g.subjind(indsort(1:min(g.methodparam, length(indsort))));
            p = histc(subjs, edges);
            if strcmpi(g.method, 'relentropy')
                p = p(1:end-1)./totcount; 
            end
            p = p/sum(p);
            p(p == 0) = [];
            for tmpi = 1:length(g.weight)
                prob3d{tmpi}(InsideIndices(1,i), InsideIndices(2,i), InsideIndices(3,i)) = -sum(p.*log(p));
            end
        else
            % Distance calculation
            ordsubjs = g.subjind(indsort);
            use_dipoles = [];
            for index = 1:length(vals)
                tmpind = find(ordsubjs == vals(index), 1);
                if ~isempty(tmpind)
                    use_dipoles(index) = tmpind;
                end
            end
            for tmpi = 1:length(g.weight)
                prob3d{tmpi}(InsideIndices(1,i), InsideIndices(2,i), InsideIndices(3,i)) = ...
                   sum(tmpweights{tmpi}(use_dipoles).*exp(-tmpsort(use_dipoles)/(2*g.methodparam^2)));
            end
        end
        
        if mod(i, 100) == 0
            fprintf('%d ', i); 
        end
    end
    fprintf('\n');
    
else 
    % 'alldistance' method - scan dipoles (more efficient)
    fprintf('Processing dipoles (of %d): ', length(allx));
    
    for tmpi = 1:length(g.weight)
        tmpprob{tmpi} = zeros(1, size(InsidePoints, 2));
    end
    
    for i = 1:length(allx)
        % Compute distances from this dipole to all inside voxels
        alldists = (InsidePoints(1,:) - allx(i)).^2 + ...
                   (InsidePoints(2,:) - ally(i)).^2 + ...
                   (InsidePoints(3,:) - allz(i)).^2;
        
        % Apply Gaussian kernel
        for tmpi = 1:length(g.weight)
            tmpprob{tmpi} = tmpprob{tmpi} + g.weight{tmpi}(i) * exp(-alldists/(2*g.methodparam^2));
            if any(isinf(tmpprob{tmpi}))
                error('Infinite value in probability calculation'); 
            end
        end
        
        if mod(i, 50) == 0
            fprintf('%d ', i); 
        end
    end
    fprintf('\n');
    
    % Copy values to 3D mesh
    fprintf('Filling 3D volume...\n');
    for i = 1:length(Inside)
        pnts = allinds(:, Inside(i));
        for tmpi = 1:length(g.weight)
            prob3d{tmpi}(pnts(1), pnts(2), pnts(3)) = tmpprob{tmpi}(i);
        end
    end
end

% ========================================================================
% STEP 9: Normalize density
% ========================================================================
if strcmpi(g.method, 'alldistance') && strcmpi(g.normalization, 'on')
    fprintf('\nNormalizing to dipoles/mm^3...\n');
    for i = 1:length(g.weight)
        if any(prob3d{i}(:) < 0)
            warning('Some probabilities are negative - normalization may be problematic');
        end
        totval = sum(prob3d{i}(:));
        totdip = length(allx);
        prob3d{i} = (prob3d{i}/totval * totdip/voxvol * 1000) / g.nsessions;
    end
    fprintf('Max density: %.4f dipoles/cc\n', max(prob3d{1}(:)));
end

% ========================================================================
% STEP 10: Resample to full resolution
% ========================================================================
if g.subsample ~= 1
    fprintf('\nResampling to full MRI resolution...\n');
    for i = 1:length(g.weight)
        prob3d{i} = prob3d{i} / g.subsample;
        newprob3d = zeros(mri.dim);
        X = ceil(mri.xgrid / g.subsample);
        Y = ceil(mri.ygrid / g.subsample);
        Z = ceil(mri.zgrid / g.subsample);
        for index = 1:size(newprob3d, 3)
            newprob3d(:,:,index) = prob3d{i}(X, Y, Z(index));
        end    
        prob3d{i} = newprob3d;
    end
end

% ========================================================================
% STEP 11: Apply 3D smoothing
% ========================================================================
if g.smooth ~= 0
    fprintf('\nApplying 3D smoothing (kernel size: %.1f)...\n', g.smooth);
    for i = 1:length(g.weight)
        prob3d{i} = smooth3d(prob3d{i}, g.smooth);
    end
end

% ========================================================================
% STEP 12: Normalize to joint probability if requested
% ========================================================================
if strcmpi(g.norm2JointProb, 'on')
    fprintf('\nNormalizing to joint probability (sum = 1)...\n');
    for i = 1:length(g.weight)
        prob3d{i} = prob3d{i} / sum(prob3d{i}(:));
    end
end

% ========================================================================
% STEP 13: Plot results
% ========================================================================
if strcmpi(g.plot, 'on')
    fprintf('\nPlotting density map...\n');
    fprintf('Max density value: %.4f\n', max(prob3d{1}(:)));
    fprintf('Min density value: %.4f\n', min(prob3d{1}(:)));
    fprintf('Mean density value: %.4f\n', mean(prob3d{1}(:)));
    
    mri3dplot(prob3d, mri, g.plotargs{:});
else
    if ishandle(gcf)
        close(gcf);
    end
end

fprintf('\nDone!\n');

end