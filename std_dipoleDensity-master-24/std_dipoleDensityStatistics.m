% dipoleDensityStatistics() - Modified to perform nonparametric
%                             dipoleDensity statistics.
%
% Usage:
%               >> pValueVoxels = dipoleDensityStatistics(dipoles1, dipoles2, 'key',val, ... );
%
%
% See also: std_dipoleDensity(), dipplot(), mri3dplot(), Fieldtrip: find_inside_vol() 
%
% Author: Arnaud Delorme & Scott Makeig, SCCN, INC, UCSD
%
% History
% 04/11/2017 Makoto. Modified for statistics.
% 02/19/2013 Makoto. 'norm2JointProb' added.

% Copyright (C) Makoto Miyakoshi, Arnaud Delorme, & Scott Makeig, SCCN/INC/UCSD
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


function [trueDifference, uncorrectedPvalues3D, correctionMask3D] = std_dipoleDensityStatistics(dipoles1, dipoles2, subjectIdx1, subjectIdx2, varargin)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Obtain optional inputs. %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g = finputcheck(varargin, { 'subjind'     'integer'  []               [];
                            'method' 'string' { 'relentropy','entropy','distance','alldistance' } 'alldistance';
                            'methodparam' 'real'     []               20; 
                            'weight'      { 'real','cell' }  []               [];
                            'smooth'      'real'     []               0;
                            'nsessions'   'integer'  []               1;
                            'subsample'   'integer'  []               2;
                            'plotargs'    'cell'     []               {};
                            'plot'        'string'  { 'on','off' }    fastif(nargout == 0, 'on', 'off');
                            'dipplot'     'string'  { 'on','off' }   'off';
                            'coordformat' 'string'  { 'mni','spherical' }   'mni';
                            'normalization' 'string'  { 'on','off' } 'on';
                            'volmesh_fname' 'string'  []  'volmesh_local.mat';
                            'mri'         { 'struct','string' } [] '';
                            'numIteration' 'real'    []        1000
                            'pValue'       'real'    []        0.05     });

                        
                        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Obtain dipole locations for the two conditions compared. %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
totalDipoleCounter = 0;
for dipoleIdx = 1:length(dipoles1)
    totalDipoleCounter = totalDipoleCounter + 1;
    allx1(totalDipoleCounter) = dipoles1(dipoleIdx).posxyz(1,1);
    ally1(totalDipoleCounter) = dipoles1(dipoleIdx).posxyz(1,2);
    allz1(totalDipoleCounter) = dipoles1(dipoleIdx).posxyz(1,3);
    alli1(totalDipoleCounter) = dipoleIdx;
    
    if size(dipoles1(dipoleIdx).posxyz, 1) == 2
        totalDipoleCounter = totalDipoleCounter + 1;
        allx1(totalDipoleCounter) = dipoles1(dipoleIdx).posxyz(2,1);
        ally1(totalDipoleCounter) = dipoles1(dipoleIdx).posxyz(2,2);
        allz1(totalDipoleCounter) = dipoles1(dipoleIdx).posxyz(2,3);
        alli1(totalDipoleCounter) = dipoleIdx;
    end
end

totalDipoleCounter = 0;
for dipoleIdx = 1:length(dipoles2)
    totalDipoleCounter = totalDipoleCounter + 1;
    allx2(totalDipoleCounter) = dipoles2(dipoleIdx).posxyz(1,1);
    ally2(totalDipoleCounter) = dipoles2(dipoleIdx).posxyz(1,2);
    allz2(totalDipoleCounter) = dipoles2(dipoleIdx).posxyz(1,3);
    alli2(totalDipoleCounter) = dipoleIdx;
    
    if size(dipoles2(dipoleIdx).posxyz, 1) == 2
        totalDipoleCounter = totalDipoleCounter + 1;
        allx2(totalDipoleCounter) = dipoles2(dipoleIdx).posxyz(2,1);
        ally2(totalDipoleCounter) = dipoles2(dipoleIdx).posxyz(2,2);
        allz2(totalDipoleCounter) = dipoles2(dipoleIdx).posxyz(2,3);
        alli2(totalDipoleCounter) = dipoleIdx;
    end
end



%%%%%%%%%%%%%%%%%%%%%%
%%% Load MRI file. %%%
%%%%%%%%%%%%%%%%%%%%%%
dipfitdefs;
load('-mat', template_models(1).mrifile); % load mri variable
g.mri = mri;

% I think the following calculation is wrong. Makoto 04/11/2017.
% % Compute voxel size.
% point1 = g.mri.transform * [ 1 1 1 1 ]';
% point2 = g.mri.transform * [ 2 2 2 1 ]';
% voxvol = sum((point1(1:3)-point2(1:3)).^2)*g.subsample^3; % in mm



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compute volume inside head mesh. %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dipfitdefs; % get the location of standard BEM volume file
tmp = load('-mat',DIPOLEDENSITY_STDBEM); % load MNI mesh
filename = g.volmesh_fname; %
disp('Loading file containing inside/outide voxel indices...');
load('-mat',filename);
InsidePoints  = allpoints(:, Inside);
InsideIndices = allinds(:,   Inside);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compute dipole density. %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tmpprob1 = zeros(1, size(InsidePoints,2));
for i = 1:size(allx1,2)
    alldists = (InsidePoints(1,:) - allx1(i)).^2 + ...
               (InsidePoints(2,:) - ally1(i)).^2 + ...
               (InsidePoints(3,:) - allz1(i)).^2;
    tmpprob1 = tmpprob1 + exp(-alldists/(2*g.methodparam^2)); % 3-D gaussian smooth, no weighting
end

tmpprob2 = zeros(1, size(InsidePoints,2));
for i = 1:size(allx2,2)
    alldists = (InsidePoints(1,:) - allx2(i)).^2 + ...
               (InsidePoints(2,:) - ally2(i)).^2 + ...
               (InsidePoints(3,:) - allz2(i)).^2;
    tmpprob2 = tmpprob2 + exp(-alldists/(2*g.methodparam^2)); % 3-D gaussian smooth, no weighting
end

trueProb3d_1 = zeros(ceil(g.mri.dim/g.subsample));
trueProb3d_2 = zeros(ceil(g.mri.dim/g.subsample));
for i = 1:length(Inside)
    pnts = allinds(:,Inside(i));
    trueProb3d_1(pnts(1), pnts(2), pnts(3)) = tmpprob1(i);
    trueProb3d_2(pnts(1), pnts(2), pnts(3)) = tmpprob2(i);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Normalize dipole density for all 4x4x4mm voxels.  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% g.mri.hdr.mat is an affin matrix, which shows 2x2x2mm.
% g.mri.hdr.dim shows 91x109x91.
% The original space coordinate must be 182x218x182.
% Now, size(prob3d_1{1}) is 46x55x46, which is a half of 91x109x91.
% Hence, 46x55x46 should correspond to 4x4x4 mm voxel size.
% The current unit is dipole density per 4x4x4 mm.
% To convert it to dipole density/cc, multiply (10^3)/(4^3).
sumDipoleDensity1 = sum(trueProb3d_1(:));  % total values in the head
sumDipoleDensity2 = sum(trueProb3d_2(:));  % total values in the head
trueProb3d_1 = trueProb3d_1*(1/sumDipoleDensity1); % Joint probability. Unit is dipole density per 4x4x4 mm voxel size
trueProb3d_2 = trueProb3d_2*(1/sumDipoleDensity2); % Joint probability.



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compute true difference. %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
trueDifference = trueProb3d_1 - trueProb3d_2;
trueDifferenceLinearized = trueDifference(:);
nonzeroVoxelIdx = find(trueDifferenceLinearized); % 30,000 x 2000 == 458 MB. Maybe acceptable.



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compute surrogate distribution. %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Identify dual dipoles
[uniqueDipoleIdx1, singleOrDualDipoles1] = unique(alli1);
[uniqueDipoleIdx2, singleOrDualDipoles2] = unique(alli2);
dualDipoleIdx1 = find(diff(singleOrDualDipoles1)==2);
dualDipoleIdx2 = find(diff(singleOrDualDipoles2)==2);

% Construct dualdipole-compatible dipole-to-subject list for dipoles1
subjectIdx1IncludingDualDipoles = zeros(length(alli1), 1);
if isempty(dualDipoleIdx1)
    % If there are no dual dipoles, use subjectIdx1 directly
    subjectIdx1IncludingDualDipoles = subjectIdx1(:);
else
    for n = 1:length(dualDipoleIdx1)
        if n == 1
            originalPart = subjectIdx1(1:dualDipoleIdx1(n));
            insertedPart = subjectIdx1(dualDipoleIdx1(n));
            subjectIdx1IncludingDualDipoles(1:dualDipoleIdx1(n)+1) = [originalPart insertedPart];
        else
            originalPart = subjectIdx1(dualDipoleIdx1(n-1)+1:dualDipoleIdx1(n));
            insertedPart = subjectIdx1(dualDipoleIdx1(n));
            subjectIdx1IncludingDualDipoles(dualDipoleIdx1(n-1)+n:dualDipoleIdx1(n)+n) = [originalPart insertedPart];
        end
    end
    % Only compute lastPiece if dualDipoleIdx1 is not empty
    if dualDipoleIdx1(end)+1 <= length(subjectIdx1)
        lastPiece = subjectIdx1(dualDipoleIdx1(end)+1:end);
        subjectIdx1IncludingDualDipoles(end-length(lastPiece)+1:end) = lastPiece;
    else
        % If dualDipoleIdx1(end)+1 exceeds length(subjectIdx1), fill remaining with last subject ID
        lastSubject = subjectIdx1(dualDipoleIdx1(end));
        subjectIdx1IncludingDualDipoles(dualDipoleIdx1(end)+length(dualDipoleIdx1)+1:end) = lastSubject;
    end
end

% Construct dualdipole-compatible dipole-to-subject list for dipoles2
subjectIdx2IncludingDualDipoles = zeros(length(alli2), 1);
if isempty(dualDipoleIdx2)
    % If there are no dual dipoles, use subjectIdx2 directly
    subjectIdx2IncludingDualDipoles = subjectIdx2(:);
else
    for n = 1:length(dualDipoleIdx2)
        if n == 1
            originalPart = subjectIdx2(1:dualDipoleIdx2(n));
            insertedPart = subjectIdx2(dualDipoleIdx2(n));
            subjectIdx2IncludingDualDipoles(1:dualDipoleIdx2(n)+1) = [originalPart insertedPart];
        else
            originalPart = subjectIdx2(dualDipoleIdx2(n-1)+1:dualDipoleIdx2(n));
            insertedPart = subjectIdx2(dualDipoleIdx2(n));
            subjectIdx2IncludingDualDipoles(dualDipoleIdx2(n-1)+n:dualDipoleIdx2(n)+n) = [originalPart insertedPart];
        end
    end
    % Only compute lastPiece if dualDipoleIdx2 is not empty
    if dualDipoleIdx2(end)+1 <= length(subjectIdx2)
        lastPiece = subjectIdx2(dualDipoleIdx2(end)+1:end);
        subjectIdx2IncludingDualDipoles(end-length(lastPiece)+1:end) = lastPiece;
    else
        % If dualDipoleIdx2(end)+1 exceeds length(subjectIdx2), fill remaining with last subject ID
        lastSubject = subjectIdx2(dualDipoleIdx2(end));
        subjectIdx2IncludingDualDipoles(dualDipoleIdx2(end)+length(dualDipoleIdx2)+1:end) = lastSubject;
    end
end
% Combine two groups.
dipole1Subj = unique(subjectIdx1);
dipole2Subj = unique(subjectIdx2);
combinedSubj = [subjectIdx1IncludingDualDipoles' subjectIdx2IncludingDualDipoles'];
combinedUniqueSubj = [dipole1Subj dipole2Subj];
combinedAllX = [allx1 allx2];
combinedAllY = [ally1 ally2];
combinedAllZ = [allz1 allz2];

% Perform permutation test.
surroDifferenceStack = zeros(length(nonzeroVoxelIdx), g.numIteration);
for iterationIdx = 1:g.numIteration
    tic
    
    % Permute data across the two conditions.
    permIdx    = randperm(length(combinedUniqueSubj));
    surroSubj1 = permIdx(1:length(dipole1Subj));
    surroSubj2 = permIdx(length(dipole1Subj)+1:end);
    
    surroIcIdx1 = find(ismember(combinedSubj, surroSubj1));
    surroIcIdx2 = find(ismember(combinedSubj, surroSubj2));
    
    surroX1 = combinedAllX(surroIcIdx1);
    surroY1 = combinedAllY(surroIcIdx1);
    surroZ1 = combinedAllZ(surroIcIdx1);
    surroX2 = combinedAllX(surroIcIdx2);
    surroY2 = combinedAllY(surroIcIdx2);
    surroZ2 = combinedAllZ(surroIcIdx2);
    
    % Compute surrogate dipole density for each condition.
    surroDipProbability1 = zeros(1, size(InsidePoints,2));
    for i = 1:size(surroX1,2)
        alldists = (InsidePoints(1,:) - surroX1(i)).^2 + ...
                   (InsidePoints(2,:) - surroY1(i)).^2 + ...
                   (InsidePoints(3,:) - surroZ1(i)).^2;
        surroDipProbability1 = surroDipProbability1 + exp(-alldists/(2*g.methodparam^2)); % 3-D gaussian smooth, no weighting
    end
    
    surroDipProbability2 = zeros(1, size(InsidePoints,2));
    for i = 1:size(surroX2,2)
        alldists = (InsidePoints(1,:) - surroX2(i)).^2 + ...
                   (InsidePoints(2,:) - surroY2(i)).^2 + ...
                   (InsidePoints(3,:) - surroZ2(i)).^2;
        surroDipProbability2 = surroDipProbability2 + exp(-alldists/(2*g.methodparam^2)); % 3-D gaussian smooth, no weighting
    end
    
    surroProb3d_1 = zeros(ceil(g.mri.dim/g.subsample));
    surroProb3d_2 = zeros(ceil(g.mri.dim/g.subsample));
    for i = 1:length(Inside)
        pnts = allinds(:,Inside(i));
        surroProb3d_1(pnts(1), pnts(2), pnts(3)) = surroDipProbability1(i);
        surroProb3d_2(pnts(1), pnts(2), pnts(3)) = surroDipProbability2(i);
    end

    % Normalize dipole densities to 4x4x4mm voxels
    surroProb3d_1 = surroProb3d_1/sum(surroProb3d_1(:));
    surroProb3d_2 = surroProb3d_2/sum(surroProb3d_2(:));
    
    % Store the result to the stack. 
    surroDifference = surroProb3d_1 - surroProb3d_2;
    surroDifferenceLinearized = surroDifference(:);
    surroDifferenceStack(:,iterationIdx) = surroDifferenceLinearized(nonzeroVoxelIdx);
    
    % Display time elapsed
    if mod(iterationIdx, 10) == 0
        timeElapsed = toc;
        timeRemaining = timeElapsed*(g.numIteration-iterationIdx);
        if     timeRemaining >= 3600
            disp(sprintf('%.0f/%.0f %.1f hours remaining.',   iterationIdx, g.numIteration, timeRemaining/3600));
        elseif timeRemaining >  60
            disp(sprintf('%.0f/%.0f %.1f minutes remaining.', iterationIdx, g.numIteration, timeRemaining/60));
        else
            disp(sprintf('%.0f/%.0f %.0f seconds remaining.', iterationIdx, g.numIteration, timeRemaining));
        end
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compute nonparametric statistics. %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
uncorrectedPvalues = stat_surrogate_pvals(surroDifferenceStack, trueDifferenceLinearized(nonzeroVoxelIdx), 'both');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compute Family-Wise Error Rate (FWER) correction. %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute critical values for Family-Wise Error Rate 
alpha = g.pValue;
maxValues = max(surroDifferenceStack,[],1);
minValues = min(surroDifferenceStack,[],1);
surroDistribution = [minValues maxValues];
criticalValues = prctile(surroDistribution, [alpha*100/2 100-(alpha*100)/2]);
correctionMask = trueDifferenceLinearized(nonzeroVoxelIdx) < criticalValues(1) | trueDifferenceLinearized(nonzeroVoxelIdx) > criticalValues(2);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Recover the 3D grid of the standard dimensions. %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
uncorrectedPvaluesLinearlized = zeros(length(trueDifferenceLinearized),1);
uncorrectedPvaluesLinearlized(nonzeroVoxelIdx) = uncorrectedPvalues;
uncorrectedPvalues3D = reshape(uncorrectedPvaluesLinearlized, size(trueDifference));
correctionMaskLinearlized = zeros(length(trueDifferenceLinearized),1);
correctionMaskLinearlized(nonzeroVoxelIdx) = double(correctionMask);
correctionMask3D = reshape(correctionMaskLinearlized, size(trueDifference));

uncorrectedPvalues3D = uncorrectedPvalues3D/g.subsample;
correctionMask3D     = correctionMask3D/g.subsample;
uncorrectedPvalues3dStandardDimensions = zeros(g.mri.dim);
correctionMask3dStandardDimensions     = zeros(g.mri.dim);
X = ceil(g.mri.xgrid/g.subsample);
Y = ceil(g.mri.ygrid/g.subsample);
Z = ceil(g.mri.zgrid/g.subsample);
for Zindex = 1:size(uncorrectedPvalues3dStandardDimensions,3)
    uncorrectedPvalues3dStandardDimensions(:,:,Zindex) = uncorrectedPvalues3D(X,Y,Z(Zindex));
    correctionMask3dStandardDimensions(:,:,Zindex)     = correctionMask3D(X,Y,Z(Zindex));
end
uncorrectedPvalues3D = uncorrectedPvalues3dStandardDimensions;
correctionMask3D     = correctionMask3dStandardDimensions;
