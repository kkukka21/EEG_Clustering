% Look up dipole structure labels from Talairach, and add to EEGLAB dataset (in the .dipfit field)
% EEG = eeg_lookup_talairach(EEG)
%
% In:
%   EEG : EEGLAB data set with .dipfit structure
%
%   ConfusionSphere : radius of assumed sphere of confusion around dipfit locations (to arrive at 
%                     probabilities), in milimeters (default: 10)
%
% Out:
%   EEG : EEGLAB data set with associated labels
%
% Example:
%   % load data set and do a lookup
%   eeg = pop_loadset('/data/projects/RSVP/exp53/realtime/exp53_target_epochs.set')
%   labeled = eeg_lookup_talairach(eeg)
%
%   % show structure labels and associated probabilities for component/dipole #17
%   labeled.dipfit.model(17).structures
%   labeled.dipfit.model(17).probabilities
%
% TODO:
%   % Replace sphere by a Gaussian with std. dev.
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2011-04-06
% History:
% 10/02/2018 Makoto. Modified. 

function [structuresFinal, probabilitiesFinal] = std_dipoleDensity_eeg_lookup_talairach(inputXyz, confusion_sphere)

if ~exist('confusion_sphere','var')
    confusion_sphere = 10; end

currentFolderPath = which('std_pop_dipoleDensity');
if ~exist('org.talairach.Database','class')
    addpathName1      = [currentFolderPath(1:end-24) filesep 'dependency' filesep 'Talairach' filesep 'talairach.jar'];
    javaaddpath(addpathName1);
end
    % if ~exist('org.talairach.Database','class')
    %     javaaddpath('/data/common/brain_atlases/Talairach/talairach.jar'); end

db = org.talairach.Database;
addpathName2 = [currentFolderPath(1:end-24) filesep 'dependency' filesep 'Talairach' filesep 'talairach.nii'];
db.load(addpathName2);

% Obtain all the labels.
p = icbm_spm2tal(inputXyz);
allLabels = cellfun(@(d)char(d),cell(db.search_range(p(1),p(2),p(3),confusion_sphere*1)),'UniformOutput',false);

% Extract Level3 and Level5 only. See http://www.talairach.org/labels.html.
level3List = cell(length(allLabels),1);
level5List = cell(length(allLabels),1);
for labelIdx = 1:length(allLabels)
    currentLabel = hlp_split(sprintf('%s,',allLabels{labelIdx,1}),',');
    level3List{labelIdx,1} = currentLabel{3};
    level5List{labelIdx,1} = currentLabel{5};
end

%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Clean the labels. %%%
%%%%%%%%%%%%%%%%%%%%%%%%%
% This is for Level3. x is to exclude.
%         Angular Gyrus
%         Anterior Cingulate
%       x Caudate
%       x Cerebellar Lingual
%       x Cerebellar Tonsil
%         Cingulate Gyrus
%       x Claustrum
%       x Culmen
%       x Culmen of Vermis
%         Cuneus
%       x Declive
%       x Declive of Vermis
%       x Extra-Nuclear
%       x Fastigium
%       x Fourth Ventricle
%         Fusiform Gyrus
%         Inferior Frontal Gyrus
%         Inferior Occipital Gyrus
%         Inferior Parietal Lobule
%       x Inferior Semi-Lunar Lobule
%         Inferior Temporal Gyrus
%         Insula
%       x Lateral Ventricle
%       x Lentiform Nucleus
%         Lingual Gyrus
%         Medial Frontal Gyrus
%         Middle Frontal Gyrus
%         Middle Occipital Gyrus
%         Middle Temporal Gyrus
%       x Nodule
%         Orbital Gyrus
%         Paracentral Lobule
%       x Parahippocampal Gyrus
%         Postcentral Gyrus
%         Posterior Cingulate
%         Precentral Gyrus
%         Precuneus
%       x Pyramis
%       x Pyramis of Vermis
%         Rectal Gyrus
%       x Subcallosal Gyrus
%       x Sub-Gyral
%         Superior Frontal Gyrus
%         Superior Occipital Gyrus
%         Superior Parietal Lobule
%         Superior Temporal Gyrus
%         Supramarginal Gyrus
%       x Thalamus
%       x Third Ventricle
%         Transverse Temporal Gyrus
%       x Tuber
%       x Tuber of Vermis
%       x Uncus
%       x Uvula
%       x Uvula of Vermis

level3ExclusionList = ...
strcmp(level3List, 'Caudate')| ...
strcmp(level3List, 'Cerebellar Lingual')| ...
strcmp(level3List, 'Cerebellar Tonsil')| ...
strcmp(level3List, 'Claustrum')| ...
strcmp(level3List, 'Culmen')| ...
strcmp(level3List, 'Culmen of Vermis')| ...
strcmp(level3List, 'Caudate')| ...
strcmp(level3List, 'Caudate')| ...
strcmp(level3List, 'Declive')| ...
strcmp(level3List, 'Declive of Vermis')| ...
strcmp(level3List, 'Extra-Nuclear')| ...
strcmp(level3List, 'Fastigium')| ...
strcmp(level3List, 'Fourth Ventricle')| ...
strcmp(level3List, 'Inferior Semi-Lunar Lobule')| ...
strcmp(level3List, 'Lateral Ventricle')| ...
strcmp(level3List, 'Lentiform Nucleus')| ...
strcmp(level3List, 'Nodule')| ...
strcmp(level3List, 'Parahippocampal Gyrus')| ...
strcmp(level3List, 'Pyramis')| ...
strcmp(level3List, 'Pyramis of Vermis')| ...
strcmp(level3List, 'Subcallosal Gyrus')| ...
strcmp(level3List, 'Sub-Gyral')| ...
strcmp(level3List, 'Thalamus')| ...
strcmp(level3List, 'Third Ventricle')| ...
strcmp(level3List, 'Tuber')| ...
strcmp(level3List, 'Tuber of Vermis')| ...
strcmp(level3List, 'Uncus')| ...
strcmp(level3List, 'Uvula')| ...
strcmp(level3List, 'Uvula of Vermis');

% Replace non-EEG anatomy labels with '*'.
[level3List{level3ExclusionList}] = deal('*');

% Extract 'Brodmann area XX' from level5List.
findBrodmannArea    = strfind(level5List, 'Brodmann area');
level5ExclusionList = cellfun('isempty', findBrodmannArea);
[level5List{level5ExclusionList}] = deal('*');

% Compute level3 probabilities.
[level3Structures, dummy, idxs] = unique(level3List);
uniqueIdxsCount = zeros(length(level3Structures),1);
for idxsIdx = 1:max(idxs)
    uniqueIdxsCount(idxsIdx) = sum(idxs==idxsIdx);
end
level3Probabilities = uniqueIdxsCount/sum(uniqueIdxsCount);
[level3Probabilities,reindex] = sort(level3Probabilities,'descend');
level3Structures = level3Structures(reindex);
mask = ~strcmp(level3Structures,'*');
level3Structures    = level3Structures(mask);
level3Probabilities = level3Probabilities(mask);

% Compute level5 probabilities.
[level5Structures, dummy, idxs] = unique(level5List);
uniqueIdxsCount = zeros(length(level5Structures),1);
for idxsIdx = 1:max(idxs)
    uniqueIdxsCount(idxsIdx) = sum(idxs==idxsIdx);
end
level5Probabilities = uniqueIdxsCount/sum(uniqueIdxsCount);
[level5Probabilities,reindex] = sort(level5Probabilities,'descend');
level5Structures = level5Structures(reindex);
mask = ~strcmp(level5Structures,'*');
level5Structures    = level5Structures(mask);
level5Probabilities = level5Probabilities(mask);

% Prepare the final output.
structuresFinal    = [level3Structures; level5Structures];
probabilitiesFinal = [level3Probabilities; level5Probabilities];



function outpoints = icbm_spm2tal(inpoints)
%
% This function converts coordinates from MNI space (normalized 
% using the SPM software package) to Talairach space using the 
% icbm2tal transform developed and validated by Jack Lancaster 
% at the Research Imaging Center in San Antonio, Texas.
%
% http://www3.interscience.wiley.com/cgi-bin/abstract/114104479/ABSTRACT
% 
% FORMAT outpoints = icbm_spm2tal(inpoints)
% Where inpoints is N by 3 or 3 by N matrix of coordinates
% (N being the number of points)
%
% ric.uthscsa.edu 3/14/07

% find which dimensions are of size 3
dimdim = find(size(inpoints) == 3);
if isempty(dimdim)
  error('input must be a N by 3 or 3 by N matrix')
end

% 3x3 matrices are ambiguous
% default to coordinates within a row
if dimdim == [1 2]
  disp('input is an ambiguous 3 by 3 matrix')
  disp('assuming coordinates are row vectors')
  dimdim = 2;
end

% transpose if necessary
if dimdim == 2
  inpoints = inpoints';
end

% Transformation matrices, different for each software package
icbm_spm = [0.9254 0.0024 -0.0118 -1.0207
	   	   -0.0048 0.9316 -0.0871 -1.7667
            0.0152 0.0883  0.8924  4.0926
            0.0000 0.0000  0.0000  1.0000];

% apply the transformation matrix
inpoints = [inpoints; ones(1, size(inpoints, 2))];
inpoints = icbm_spm * inpoints;

% format the outpoints, transpose if necessary
outpoints = inpoints(1:3, :);
if dimdim == 2
  outpoints = outpoints';
end



function res = hlp_split(str,delims)
% Split a string according to some delimiter(s).
% Result = hlp_split(String,Delimiters)
%
% In:
%   String : a string (char vector)
%
%   Delimiters : a vector of delimiter characters (includes no special support for escape sequences)
%
% Out:
%   Result : a cell array of (non-empty) non-Delimiter substrings in String
%
% Examples:
%   % split a string at colons and semicolons; returns a cell array of four parts
%   hlp_split('sdfdf:sdfsdf;sfdsf;;:sdfsdf:',':;')
% 
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-11-05

pos = find(diff([0 ~sum(bsxfun(@eq,str(:)',delims(:)),1) 0]));
res = cell(~isempty(pos),length(pos)/2);
for k=1:length(res)
    res{k} = str(pos(k*2-1):pos(k*2)-1); end
