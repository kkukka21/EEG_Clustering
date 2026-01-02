% std_dipoleDensity(): plots dipoles and dipole density. Requies
%                      inputs from std_pop_dipplotWithDensity. In the
%                      output, the unit measure is joint probability
%                      (i.e. sum of all voxel values == 1).
% Usage:    
%                >> std_dipoleDensity(STUDY, ALLEEG, varargin);  
% Inputs:
%   STUDY    - EEGLAB STUDY set comprising some or all of the EEG datasets in ALLEEG.
%   ALLEEG   - global EEGLAB vector of EEG structures for the dataset(s) included in 
%              the STUDY. ALLEEG for a STUDY set is typically created using load_ALLEEG().  
%   varargin - user input from std_pop_dipplotWithDensity.
%
% See also: eegplugin_std_dipoleDensity() std_pop_dipplotWithDensity() dipplot() mri3dplot()        

% Author: Makoto Miyakoshi, JSPS/SCCN, INC, UCSD. Cincinnati Children's Hospital
%         Luca Pion-Tonacini, SCCN, INC, UCSD.
%         compute_centroid written by Hilit Serby, Arnaud Delorme, Scott Makeig, SCCN, INC, UCSD, June, 2005
%
% History:
% 09/18/2024 Makoto. When STUDY.design.variable is empty, generate dummy variables for group labels.
% 08/23/2024 Makoto and Komal. Re-visiting this plugin to make it work again.
% 10/02/2018 Makoto. Dipole probability density calculation fixed.
% 05/22/2017 Makoto. Fixed the calculation of a cluster centroid when bilateral dipoles are present.
% 03/16/2017 Makoto. Fixed the bilateral dipole case to choose the consistent one with cluster centroid.
% 03/06/2017 Makoto. Changed to FWHM.
% 09/13/2016 ver 0.25 by Makoto. Standard deviation and error for dipole centroid supported. 
% 10/29/2015 ver 0.24 by Makoto. 'session' supported.
% 05/20/2015 ver 0.23 by Makoto. private function 'hlp_varargin2struct' added (Thanks Jens Bernhardsson)
% 01/27/2015 ver 0.21 by Makoto. 'group' could be non-present, var1, or var2.
% 01/16/2015 ver 0.20 by Luca and Makoto. brainBlobBrowser added.
% 05/22/2013 ver 3.3 by Makoto. lightangle for sagittal and colonal views.
% 04/01/2013 ver 3.2 by Makoto. Color scheme default optimized.
% 03/29/2013 ver 3.1 by Makoto. Added cmin cmax. mir3dplot() is fixed accordingly.
% 03/20/2013 ver 3.0 by Makoto. Talairach coordinate of the blob peak (difference plot only) output. Color scheme improved.
% 03/11/2013 ver 2.3 by Makoto. Bug fixed (specifying ranges in the difference plot caused scale) 
% 02/19/2013 ver 2.2 by Makoto. norm2JointProb added.
% 02/06/2013 ver 2.1 by Makoto. Color scale upper limit added.
% 01/24/2013 ver 2.0 by Makoto. Difference between groups supported.
% 12/03/2012 ver 1.1 by Makoto. Save figure added.
% 11/23/2012 ver 1.0 by Makoto. Created.

% Copyright (C) 2012, Makoto Miyakoshi JSPS/SCCN,INC,UCSD
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

function std_dipoleDensity(STUDY, ALLEEG, varargin)
varargin = varargin{1,1};

currentDesign = STUDY.currentdesign;

allSource       = struct([]);
allDipName      = [];
allDipColor     = [];
allDipPlotTitle = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% find out if 'group' or 'session' is non-exist, var1, or var2 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
designLabel = {STUDY.design(currentDesign).variable.label};
groupFieldIdx = find(strcmp(designLabel, 'group')|strcmp(designLabel, 'session')); % find the slot index that has 'group' or 'session'
if isempty(groupFieldIdx) % no 'group'
    singletonFieldIdx = find(size(STUDY.cluster(1,2).setinds)==1); % find the slot that is NOT the within-subject condition
    if length(singletonFieldIdx)>1 % this means both var1 and var2 are empty
        singletonFieldIdx = 1;
    end
    groupFieldIdx = singletonFieldIdx;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% loop for all 5 selections %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for dipDensityPlotIdx = 1:5
    if varargin{1,dipDensityPlotIdx}.cluster>1 & varargin{1,dipDensityPlotIdx}.group>1
        currentCluster   = varargin{1,dipDensityPlotIdx}.cluster-1;
        currentGroup     = varargin{1,dipDensityPlotIdx}.group;
        currentColor     = varargin{1,dipDensityPlotIdx}.color;
        currentColorName = varargin{1,dipDensityPlotIdx}.colorName;
        
        %%%%%%%%%%%%%%%%%%%%%%%
        %%% separate groups %%%
        %%%%%%%%%%%%%%%%%%%%%%%
        STUDY = std_groupDipSeparator(STUDY, ALLEEG, currentCluster);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% select specific/all group(s) %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if currentGroup == 2 % show all groups
            tmpGroupName = 'all';
            source  = struct([]);
            dipName = [];
            for groupIdx = 1:length(STUDY.design(currentDesign).variable(1,groupFieldIdx).value);
                source   = [source STUDY.cluster(1,currentCluster).groupDipModels{1,groupIdx}];
                dipName  = [dipName STUDY.cluster(1,currentCluster).groupDipNames{1,groupIdx}];
            end
        else % selecting a group
            tmpGroupName = STUDY.design(currentDesign).variable(1,groupFieldIdx).value{1,currentGroup-2};
            % tmpGroupName = [tmpGroupName;repmat({' and '},1,size(tmpGroupName,2))];
            % tmpGroupName = [tmpGroupName{1:end-1}];
            source  = STUDY.cluster(1,currentCluster).groupDipModels{1,currentGroup-2};
            dipName = STUDY.cluster(1,currentCluster).groupDipNames{1,currentGroup-2};
        end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Add centroid in the end %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    
    % Check if there is dual dipoles.
    coordinates  = {source.posxyz}';
    isSingleDipole = zeros(length(coordinates),1);
    for icIdx = 1:length(coordinates)
        isSingleDipole(icIdx) = size(coordinates{icIdx},1)==1;
    end
    
    % Compute the centoid only using single dipoles. 
    centroid = computecentroid(source(find(isSingleDipole)));

    % Re-compute the centroid if dual dipole is present by using one of the dual dipoles on the side the cluster centroid exists (discussed with Scott). 
    if any(isSingleDipole == 0)
        singleDipoleSource = source;
        dualDipoleIdx = find(isSingleDipole == 0);
        for dualDipoleIdxIdx = 1:length(dualDipoleIdx)
            xCoordinates = singleDipoleSource(dualDipoleIdx(dualDipoleIdxIdx)).posxyz(:,1);
            if     centroid.posxyz(1) > 0 % The cluster centroid is in the right hemisphere.
                   singleDipoleSource(dualDipoleIdx(dualDipoleIdxIdx)).posxyz = singleDipoleSource(dualDipoleIdx(dualDipoleIdxIdx)).posxyz(find(xCoordinates>0),:);
                   singleDipoleSource(dualDipoleIdx(dualDipoleIdxIdx)).momxyz = singleDipoleSource(dualDipoleIdx(dualDipoleIdxIdx)).momxyz(find(xCoordinates>0),:);
            elseif centroid.posxyz(1) < 0 % The cluster centroid is in the left hemisphere.
                   singleDipoleSource(dualDipoleIdx(dualDipoleIdxIdx)).posxyz = singleDipoleSource(dualDipoleIdx(dualDipoleIdxIdx)).posxyz(find(xCoordinates<0),:);
                   singleDipoleSource(dualDipoleIdx(dualDipoleIdxIdx)).momxyz = singleDipoleSource(dualDipoleIdx(dualDipoleIdxIdx)).momxyz(find(xCoordinates<0),:);
            end
        end
        
        % Re-compute the centroid using one of the bilateral dipoles.
        centroid = computecentroid(singleDipoleSource);
        warning(sprintf('\nICs with (symmetric) bilateral dipoles detected.\nThe cluster centroid is computed by using only one of the bilateral dipoles on the side the cluster centroid exists.'))
    else
        singleDipoleSource = source;
    end
    source(end + 1) = centroid;
    centroidTal = mni2tal(centroid.posxyz);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Compute standard deviation and standard error of dipole centroid %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    coordinateTable = cell2mat({singleDipoleSource.posxyz}');
    coordinateTable = coordinateTable(1:end-1,:); % Exclude centroid
    standardDeviation = std(coordinateTable);
    standardError     = std(coordinateTable)/sqrt(size(coordinateTable,1));
    clusterReport     = sprintf('\nCluster: %.0f\nCentroid in MNI:    [%2.0f %2.0f %2.0f]\nStandard Deviation: [%2.0f %2.0f %2.0f]\nStandard Error    : [%2.0f %2.0f %2.0f]\n',...
        currentCluster,...
        centroid.posxyz(1),   centroid.posxyz(2),   centroid.posxyz(3), ...
        standardDeviation(1), standardDeviation(2), standardDeviation(3),...
        standardError(1),     standardError(2),     standardError(3));

    %%%%%%%%%%%%%%%%%%%%%%
    %%% prepare colors %%%
    %%%%%%%%%%%%%%%%%%%%%%
    dipColor = cell(1, length(source));
    dipColor(1:length(source)-1) = {currentColor};
    dipColor(end) = {[1 0 0]};
    
    %%%%%%%%%%%%%%%%%%%%%
    %%% prepare names %%%
    %%%%%%%%%%%%%%%%%%%%%
    centroidName = [STUDY.cluster(1,currentCluster).name ' mean'];
    dipName = [dipName centroidName];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% plot dipoledensity %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    customJet = jet(128);
    customJet(1,:)=0.3; % this is to make the background black
    if isempty(varargin{1,8})
        cmin = 0;
        cmax = [];
    else
        cmin = varargin{1,8}(1);
        cmax = varargin{1,8}(2);
    end
    
    if     varargin{1,6} == 1; % axial
        plotargs = {'mriview', 'top',  'mrislices', -30:10:70, 'cmap', customJet, 'mixfact', 0.65, 'cmin', cmin, 'cmax', cmax};
    elseif varargin{1,6} == 2; % sagittal
        plotargs = {'mriview', 'side', 'mrislices', -70:10:70, 'cmap', customJet, 'mixfact', 0.65, 'cmin', cmin, 'cmax', cmax};
    elseif varargin{1,6} == 3; % coronal
        plotargs = {'mriview', 'rear', 'mrislices', -90:10:60, 'cmap', customJet, 'mixfact', 0.65, 'cmin', cmin, 'cmax', cmax};
    end
    
    [dens3d mri] = dipoledensity(source, 'coordformat', ALLEEG(1,1).dipfit.coordformat, ...
        'methodparam', varargin{1,7}, 'plot', 'on', 'norm2JointProb', 'on', 'plotargs', plotargs);
    
    % Clarified the ambiguity of unit. 10/02/2018 Makoto.
    disp(sprintf('CORRECTION: Brightest color denotes a dipole probability density per voxel of: %.7f (Total 249728 2x2x2 mm inside-brain voxels)', max(dens3d{1}(:))))
    
    
    
    h1 = gcf;
    if iscell(tmpGroupName)
        tmpCell = tmpGroupName;
        tmpCell(2,:) = {' & '};
        tmpCell{2,end} = '';
        tmpString = [tmpCell{:}];
        tmpGroupName = tmpString;
    end
    set(h1, 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off', 'menu', 'none', 'NumberTitle','off','Name', ['Cls ' num2str(currentCluster) ' Group ' tmpGroupName  '; std_dipoleDensity()']);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Call BrainBlobBrowser %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    brainBlobBrowser('data', dens3d{1}, 'clusterReport', clusterReport, 'clusterStd', standardDeviation);
    
    allGroupDens3d{currentGroup} = dens3d;
    h1 = gcf;
    set(h1, 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off', 'menu', 'none', 'NumberTitle','off','Name', ['Cls ' num2str(currentCluster) ' Group ' tmpGroupName  '; BrainBlobBrowser()']);
    
    allSource   = [allSource source];
    allDipName  = [allDipName dipName];
    allDipColor = [allDipColor dipColor];
    allDipPlotTitle = [allDipPlotTitle currentColorName ', Cls ' num2str(currentCluster) ' Group ' tmpGroupName '; '];
    
    if varargin{1,12} == 1
        print(h1, '-deps', ['dipDensity_' allDipPlotTitle(1:end-2)], '-loose')
    end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% plot group difference with permutation test %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if varargin{9}>1 && varargin{10}>1 && ~isempty(varargin{11})
    
    % Compute surrogate difference.
    n = 1;
    currentCluster  = varargin{1,n}.cluster-1;
    currentGroup    = varargin{1,n}.group;
    dipoles1    = STUDY.cluster(1,currentCluster).groupDipModels{1,currentGroup-2};
    subjectIdx1 = STUDY.cluster(1,currentCluster).setinds{1,currentGroup-2};
    
    n = 2;
    currentCluster  = varargin{1,n}.cluster-1;
    currentGroup    = varargin{1,n}.group;    
    dipoles2    = STUDY.cluster(1,currentCluster).groupDipModels{1,currentGroup-2};
    subjectIdx2 = STUDY.cluster(1,currentCluster).setinds{1,currentGroup-2};
    
    num2str(varargin{11})
    
    [trueDifference, uncorrectedPvalues3D, correctionMask3D] = std_dipoleDensityStatistics(dipoles1, dipoles2, subjectIdx1, subjectIdx2, 'coordformat', ALLEEG(1,1).dipfit.coordformat, 'methodparam', varargin{1,7}, 'pValue', varargin{11}/100, 'numIteration', 500);

    % Mask zero values including outside the brain.
    trueDifference(trueDifference==0) = NaN;
    uncorrectedpvalueMask3D = uncorrectedPvalues3D < varargin{11}/100;
    uncorrectedPvalues3D = uncorrectedPvalues3D.*uncorrectedpvalueMask3D;
    uncorrectedPvalues3D(uncorrectedPvalues3D==0) = NaN;
    correctionMask3D(correctionMask3D==0)         = NaN;
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Plot true difference map %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cmax = 2*max(abs(trueDifference(:)));
    yTickLabel = linspace(-max(abs(trueDifference(:))), max(abs(trueDifference(:))), 5);
    mri3dplot(trueDifference, mri, plotargs{:}, 'cmax', cmax, 'mixfact', 0.45)
    set(gca, 'YTickLabel', yTickLabel)
    set(get(gca, 'title'), 'string', 'p-values (uncorr.)', 'color', [0.7 0.7 0.7])
    set(gcf, 'NumberTitle', 'off', 'Name', 'p-values (uncorr.)', 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off', 'menu', 'none');

    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Plot uncorrected p-value map %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    logUncorrPval3D = abs(log10(uncorrectedPvalues3D));
    cmax = max(logUncorrPval3D(:));
    yTickLabel = round(logspace(log10(max(nonzeros(uncorrectedPvalues3D(:)))), log10(min(nonzeros(uncorrectedPvalues3D(:)))), 5)*10000)/10000;
    mri3dplot(logUncorrPval3D, mri, plotargs{:}, 'cmax', cmax, 'mixfact', 0.45)
    set(gca, 'YTickLabel', yTickLabel)
    set(get(gca, 'title'), 'string', 'p-values (uncorr.)', 'color', [0.7 0.7 0.7])
    set(gcf, 'NumberTitle', 'off', 'Name', 'p-values (uncorr.)', 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off', 'menu', 'none');

    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Plot corrected p-value map %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    logCorrPval3D = abs(log10(uncorrectedPvalues3D.*correctionMask3D));
    mri3dplot(logCorrPval3D, mri, plotargs{:}, 'cmax', cmax, 'mixfact', 0.45)
    set(gca, 'YTickLabel', yTickLabel)
    set(get(gca, 'title'), 'string', 'p-values (FWER corr.)', 'color', [0.7 0.7 0.7])
    set(gcf, 'NumberTitle', 'off', 'Name', 'p-values (FWER corr.)', 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off', 'menu', 'none')
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Find cluster peaks and their coordinates %%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    logUncorrPval3DBinary = ~isnan(logUncorrPval3D);
    [segmentedBlob, numBlobs] = bwlabeln(logUncorrPval3DBinary);
    
    blobStats = zeros(numBlobs,2);
    for n = 1:numBlobs
        blobStats(n,:) = [n sum(segmentedBlob(:)==n)*(2^3)/1000]; % converting 1x1x1 mm to 2x2x2 mm, then mm^3 to cm^3
    end
    sortedBlobSizes = sort(blobStats(:,2), 'descend');
    
    
    
    
    for n = 1:blobIdx
        linearInd   = find(segmentedBlob==n);
        tmpBlobSize = length(linearInd);
        tmpBlob     = tmpDiffMaskRawLin(linearInd);
        % Is this blob positive or negative?
        if max(tmpBlob) > rightCut
            [peakVal,idx] = max(tmpBlob);
            blobSign      = 'Positive';
        else
            [peakVal,idx] = min(tmpBlob);
            blobSign      = 'Negative';
        end
        % store peaks of blobs
        [peakVoxel(1,1), peakVoxel(1,2), peakVoxel(1,3)] = ind2sub(size(tmpDiffMaskRaw), linearInd(idx));
        % convert values into MNI coordinate
        peakMNI = peakVoxel*mri.transform([1:3],[1:3]) + mri.transform([1:3],4)';
        % convert MNI into Talairach coordinate (DIPFIT default)
        peakTal = round(mni2tal(peakMNI));
        % display results
        disp(['Blob' num2str(n) '(' blobSign '): Size ' num2str(tmpBlobSize/1000) 'cc, peak at [' num2str(peakTal) '] (Talairach)' ])
    end
    
    
    
    tmpDiffMaskRaw     = tmpDiffMasked{1,1};
    tmpDiffMaskRawLin  = tmpDiffMaskRaw(:);
    tmpDiffMaskBinary  = ~isnan(tmpDiffMaskRaw);
    [segmentedBlob,blobIdx] = bwlabeln(tmpDiffMaskBinary);
    for n = 1:blobIdx
        linearInd   = find(segmentedBlob==n);
        tmpBlobSize = length(linearInd);
        tmpBlob     = tmpDiffMaskRawLin(linearInd);
        % Is this blob positive or negative?
        if max(tmpBlob) > rightCut
            [peakVal,idx] = max(tmpBlob);
            blobSign      = 'Positive';
        else
            [peakVal,idx] = min(tmpBlob);
            blobSign      = 'Negative';
        end
        % store peaks of blobs
        [peakVoxel(1,1), peakVoxel(1,2), peakVoxel(1,3)] = ind2sub(size(tmpDiffMaskRaw), linearInd(idx));
        % convert values into MNI coordinate
        peakMNI = peakVoxel*mri.transform([1:3],[1:3]) + mri.transform([1:3],4)';
        % convert MNI into Talairach coordinate (DIPFIT default)
        peakTal = round(mni2tal(peakMNI));
        % display results
        disp(['Blob' num2str(n) '(' blobSign '): Size ' num2str(tmpBlobSize/1000) 'cc, peak at [' num2str(peakTal) '] (Talairach)' ])
    end
    
    
        %     tmpDiffMaskRaw     = tmpDiffMasked{1,1};
        %     tmpDiffMaskRawLin  = tmpDiffMaskRaw(:);
        %     tmpDiffMaskBinary  = ~isnan(tmpDiffMaskRaw);
        %     [segBlob,blobNum] = bwlabeln(tmpDiffMaskBinary);
        %     for n = 1:blobNum
        %         linearInd   = find(segBlob==n);
        %         tmpBlobSize = length(linearInd);
        %         tmpBlob     = tmpDiffMaskRawLin(linearInd);
        %         % Is this blob positive or negative?
        %         if max(tmpBlob) > rightCut
        %             [peakVal,idx] = max(tmpBlob);
        %             blobSign      = 'Positive';
        %         else
        %             [peakVal,idx] = min(tmpBlob);
        %             blobSign      = 'Negative';
        %         end
        %         % store peaks of blobs
        %         [peakVoxel(1,1), peakVoxel(1,2), peakVoxel(1,3)] = ind2sub(size(tmpDiffMaskRaw), linearInd(idx));
        %         % convert values into MNI coordinate
        %         peakMNI = peakVoxel*mri.transform([1:3],[1:3]) + mri.transform([1:3],4)';
        %         % convert MNI into Talairach coordinate (DIPFIT default)
        %         peakTal = round(mni2tal(peakMNI));
        %         % display results
        %         disp(['Blob' num2str(n) '(' blobSign '): Size ' num2str(tmpBlobSize/1000) 'cc, peak at [' num2str(peakTal) '] (Talairach)' ])
        %     end
    

    
    % plot
    mri3dplot(tmpDiffMasked, mri, plotargs{:}, 'cmax', cmax, 'mixfact', 0.45) % this is using cmax twice but it seems ok
    disp('Oops, the previous ''Brightest color denotes a density of...'' was wrong.')
    disp('See the color bar in the figure for the correct scale.')
    set(gca, 'YTickLabel',yTickLabel)
    set(gcf, 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off', 'menu', 'none', 'Name', ['Thresholded at ' num2str(varargin{11}) '%; std_dipoleDensity()'], 'NumberTitle','off')
    
    % print out (option)
    if varargin{1,12} == 1
        print(gcf, '-deps', ['dipDensDiff_' allDipPlotTitle(1:end-2)], '-loose')
    end
end


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%% plot group difference %%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if varargin{9}>1 && varargin{10}>1 && ~isempty(varargin{11})
%     % perform subtraction
%     tmpDiff = allGroupDens3d{1,varargin{9}}{1,1}-allGroupDens3d{1,varargin{10}}{1,1};
% 
%     % mask outside brain
%     tmpDiff(tmpDiff==0)=NaN;
%     
%     % color bar range
%     if abs(max(tmpDiff(:)))>abs(min(tmpDiff(:)))
%         cmax = 2*abs(max(tmpDiff(:)));
%         yTickLabel = linspace(-abs(max(tmpDiff(:))), abs(max(tmpDiff(:))), 5);
%         yTickLabel = round(yTickLabel*1000000)/1000000;
%     else
%         cmax = 2*abs(min(min(min(tmpDiff))));
%         yTickLabel = linspace(-abs(min(tmpDiff(:))), abs(min(tmpDiff(:))), 5);
%         yTickLabel = round(yTickLabel*1000000)/1000000;
%     end
% 
%     % threshold
%     tmpDiffSort = sort(tmpDiff(:));
%     tmpDiffSort(isnan(tmpDiffSort))=[];
%     leftCut  = tmpDiffSort(round(length(tmpDiffSort)*varargin{11}/2/100));
%     rightCut = tmpDiffSort(round(length(tmpDiffSort)*(1-varargin{11}/2/100)));
%     leftSide = tmpDiff;
%     leftSide(leftSide>leftCut) = 0;
%     rightSide = tmpDiff;
%     rightSide(rightSide<rightCut) = 0;
%     tmpDiffMasked = leftSide + rightSide;
%     tmpDiffMasked(tmpDiffMasked==0) = NaN;
%     minColorResolution = (abs(max(tmpDiffMasked(:))) + abs(min(tmpDiffMasked(:))))/127; % this is to avoid 0 that is for non-significant regions
%     tmpDiffMasked = tmpDiffMasked - min(tmpDiffMasked(:)) + minColorResolution;
%     tmpDiffMasked = {tmpDiffMasked};
% 
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %%% find blobs and their peaks in the Talairach coordinate %%%
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     tmpDiffMaskRaw     = tmpDiffMasked{1,1};
%     tmpDiffMaskRawLin  = tmpDiffMaskRaw(:);
%     tmpDiffMaskBinary  = ~isnan(tmpDiffMaskRaw);
%     [segBlob,blobNum] = bwlabeln(tmpDiffMaskBinary);
%     for n = 1:blobNum
%         linearInd   = find(segBlob==n);
%         tmpBlobSize = length(linearInd);
%         tmpBlob     = tmpDiffMaskRawLin(linearInd);
%         % Is this blob positive or negative?
%         if max(tmpBlob) > rightCut
%             [peakVal,idx] = max(tmpBlob);
%             blobSign      = 'Positive';
%         else
%             [peakVal,idx] = min(tmpBlob);
%             blobSign      = 'Negative';
%         end
%         % store peaks of blobs
%         [peakVoxel(1,1), peakVoxel(1,2), peakVoxel(1,3)] = ind2sub(size(tmpDiffMaskRaw), linearInd(idx));
%         % convert values into MNI coordinate
%         peakMNI = peakVoxel*mri.transform([1:3],[1:3]) + mri.transform([1:3],4)';
%         % convert MNI into Talairach coordinate (DIPFIT default)
%         peakTal = round(mni2tal(peakMNI));
%         % display results
%         disp(['Blob' num2str(n) '(' blobSign '): Size ' num2str(tmpBlobSize/1000) 'cc, peak at [' num2str(peakTal) '] (Talairach)' ])
%     end
%     
% 
%     
%     % plot
%     mri3dplot(tmpDiffMasked, mri, plotargs{:}, 'cmax', cmax, 'mixfact', 0.45) % this is using cmax twice but it seems ok
%     disp('Oops, the previous ''Brightest color denotes a density of...'' was wrong.')
%     disp('See the color bar in the figure for the correct scale.')
%     set(gca, 'YTickLabel',yTickLabel)
%     set(gcf, 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off', 'menu', 'none', 'Name', ['Thresholded at ' num2str(varargin{11}) '%; std_dipoleDensity()'], 'NumberTitle','off')
%     
%     % print out (option)
%     if varargin{1,12} == 1
%         print(gcf, '-dpsc2', ['dipDensDiff_' allDipPlotTitle(1:end-2)], '-loose')
%     end
% end

%%%%%%%%%%%%%%%%%%%%
%%% plot dipoles %%%
%%%%%%%%%%%%%%%%%%%%
if     varargin{1,6} == 1; % axial
    viewAngle = [0 0 1];
elseif varargin{1,6} == 2; % sagittal
    viewAngle = [1 0 0];
elseif varargin{1,6} == 3; % coronal
    viewAngle = [0 -1 0];
end

dipplot(allSource, 'mri', ALLEEG(1,1).dipfit.mrifile, ...
              'meshdata', ALLEEG(1,1).dipfit.hdmfile, ...
           'coordformat', ALLEEG(1,1).dipfit.coordformat, ...
          'dipolelength', 0, 'spheres', 'on', 'projlines', 'off',...
          'view', viewAngle, 'projimg', 'off', 'color', allDipColor,...
          'dipnames', allDipName);
h2 = gcf;
set(h2, 'PaperPositionMode', 'auto', 'InvertHardcopy', 'off', 'menu', 'none', 'NumberTitle','off','Name', [allDipPlotTitle 'std_myDipPlot()']);

% adjust lightings 
if     varargin{1,6} == 2; % sagittal
    lightangle(45,180);
elseif varargin{1,6} == 3; % coronal
    lightangle(60,180);
end

if varargin{1,12} == 1
    print(h2, '-deps', ['dip_' allDipPlotTitle(1:end-2)], '-loose')
end


function dipole = computecentroid(alldipoles)
max_r = 0;
len = length(alldipoles);
dipole.posxyz = [ 0 0 0 ];
dipole.momxyz = [ 0 0 0 ];
dipole.rv = 0;
ndip = 0;
count = 0;
numNaN = 0;
warningon = 1;
for k = 1:len 
    if size(alldipoles(k).posxyz,1) == 2
        if all(alldipoles(k).posxyz(2,:) == [ 0 0 0 ])
            alldipoles(k).posxyz(2,:) = [];
            alldipoles(k).momxyz(2,:) = [];
        end;
    end;
    if ~isempty(alldipoles(k).posxyz)
        dipole.posxyz = dipole.posxyz + mean(alldipoles(k).posxyz,1);
        dipole.momxyz = dipole.momxyz + mean(alldipoles(k).momxyz,1);
        if ~isnan(alldipoles(k).rv)
            dipole.rv = dipole.rv + alldipoles(k).rv;
        else
            numNaN = numNaN+1;
        end
        count = count+1;
    elseif warningon
        disp('Some components do not have dipole information');
        warningon = 0;
    end;
end
dipole.posxyz = dipole.posxyz/count;
dipole.momxyz = dipole.momxyz/count;
dipole.rv     = dipole.rv/(count-numNaN);
if isfield(alldipoles, 'maxr')
    dipole.maxr = alldipoles(1).max_r;
end;
        
    
    
function STUDY = std_groupDipSeparator(STUDY, ALLEEG, clusterIndex)

% Reserve the space for the output.
STUDY.cluster(1,clusterIndex).groupDipModels = struct([]);
groupDipModel                                = struct([]);

% Identify where 'group' labels are stored. I write this part based on prediction of current STUDY's behavior. (08/23/2024 Makoto)
variableLabels  = cellfun(@(x) x, {STUDY.design(STUDY.currentdesign).variable.label}, 'UniformOutput', false);
groupLabelIndex = contains(variableLabels, 'group');

% Obtain all 'group' labels.
groupLabelCells = cellfun(@(x) x{1}, STUDY.design(STUDY.currentdesign).variable(groupLabelIndex).value, 'UniformOutput', false);

% If no group name specified, auto-generate dummy variables.
if strcmp(groupLabelCells, '(autoGeneratedDummy)')
    for eegIdx = 1:length(ALLEEG)
        ALLEEG(1,eegIdx).group = {'(autoGeneratedDummy)'};
    end
end

% Obtain ALLEEG's group indices.
alleegGroupCells = cellfun(@(x) x{1}, {ALLEEG.group}, 'UniformOutput', false)';

for groupIdx = 1:length(groupLabelCells)

    % Obtain the target group label.
    currentGroupLabel = groupLabelCells(groupIdx);

    % Determine who belong to the target label.
    currentGroupSetIdx = find(contains(alleegGroupCells, currentGroupLabel));
    
    % Loop for all subjects witin the group.
    groupDipModel = [];
    groupDipName  = [];
    setIcIdx      = 0;
    for setIdx = 1:length(currentGroupSetIdx)

        currentSetIdx = currentGroupSetIdx(setIdx);
        targetSetsIdxVec = find(STUDY.cluster(clusterIndex).sets==currentSetIdx);
        targetIcIdxVec   = STUDY.cluster(clusterIndex).comps(targetSetsIdxVec);

        % Loop for all ICs within the subject.
        for icIdx = 1:length(targetIcIdxVec)

            currentIcIdx  = targetIcIdxVec(icIdx);
            tmpDipModel   = ALLEEG(1,currentSetIdx).dipfit.model(1,currentIcIdx);
            tmpDipName    = sprintf('%s, IC%d', ALLEEG(1,currentSetIdx).subject, currentIcIdx);

            setIcIdx = setIcIdx+1;
            groupDipModel(1, setIcIdx).posxyz = tmpDipModel.posxyz;
            groupDipModel(1, setIcIdx).momxyz = tmpDipModel.momxyz;
            groupDipModel(1, setIcIdx).rv     = tmpDipModel.rv;
            groupDipName{ 1, setIcIdx}        = tmpDipName;
        end
    end
    STUDY.cluster(1,clusterIndex).groupDipModels{1,groupIdx} = groupDipModel;
    STUDY.cluster(1,clusterIndex).groupDipNames{ 1,groupIdx} = groupDipName;
end

% designLabel = {STUDY.design(STUDY.currentdesign).variable.label};
% groupFieldIdx = find(strcmp(designLabel, 'group')|strcmp(designLabel, 'session')); % is 'group' or 'session' non-exist, var1, or var2?
% if isempty(groupFieldIdx) % no 'group'
%     singletonFieldIdx = find(size(STUDY.cluster(1,clusterIndex).setinds)==1); % find the slot that is NOT the within-subject condition
%     if length(singletonFieldIdx)>1 % this means both var1 and var2 are empty
%         singletonFieldIdx = 1;
%     end
%     groupFieldIdx = singletonFieldIdx;
% end
% 
% currentDesign = STUDY.currentdesign;
% for groupIndex = 1:size(STUDY.cluster(1,clusterIndex).setinds, groupFieldIdx)
%     if groupFieldIdx == 1 % var1 == 'group'
%         tmpSetInds = STUDY.cluster(1,clusterIndex).setinds(groupIndex,1);
%     else                  % var2 == 'group'
%         tmpSetInds = STUDY.cluster(1,clusterIndex).setinds(1,groupIndex);
%     end
%     tmpSetInds = tmpSetInds{1,1};
% 
%     for nthTmpSetInds = 1:length(tmpSetInds)
%         trueSetIndex  = STUDY.design(currentDesign).cell(tmpSetInds(nthTmpSetInds)).dataset;
%         if groupFieldIdx == 1 % var1 == 'group'
%             tmpIcIndex = STUDY.cluster(1,clusterIndex).allinds{groupIndex,1}(1, nthTmpSetInds);
%         else % var2 == 'group'
%             tmpIcIndex = STUDY.cluster(1,clusterIndex).allinds{1,groupIndex}(1, nthTmpSetInds);
%         end
%         tmpDipModel   = ALLEEG(1,trueSetIndex).dipfit.model(1,tmpIcIndex);
%         tmpDipName    = [ALLEEG(1,trueSetIndex).subject ', IC' num2str(tmpIcIndex)];
% 
%         groupDipModel(1, nthTmpSetInds).posxyz = tmpDipModel.posxyz;
%         groupDipModel(1, nthTmpSetInds).momxyz = tmpDipModel.momxyz;
%         groupDipModel(1, nthTmpSetInds).rv     = tmpDipModel.rv;
%         groupDipName{1,  nthTmpSetInds}        = tmpDipName;
%     end
%     STUDY.cluster(1,clusterIndex).groupDipModels{1,groupIndex} = groupDipModel;
%     STUDY.cluster(1,clusterIndex).groupDipNames{1,groupIndex} = groupDipName;
%     clear groupDipModel groupDipName
% end