% plotDipoles() - A simple multiple dipole plotting tool. Basically a
%                 wrapper for dipplot().
%
% Input:
%        xyzMatrix: [nx3] Takes XYZ coordinates in MNI.
%        rgbMatrix: [nx3] Takes RGB triplets.
%        mriData:   (Optional) MRI data structure for background overlay (default: none).
%
% Output:
%       dipPlotHandle: Handle to the axes object containing the dipole plot.
%
% History:
% 05/03/2024 Makoto. Created upon the request of Komal.
% [Date] Grok. Updated to add MRI support and remove internal figure creation.

% Copyright (C) 2024 Makoto Miyakoshi. Cincinnati Children's Hospital Medical Center.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
%
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright notice,
% this list of conditions and the following disclaimer in the documentation
% and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF
% THE POSSIBILITY OF SUCH DAMAGE.

function dipPlotHandle = plotDipoles(xyzMatrix, rgbMatrix, mriData)
    % Input validation
    if nargin < 2
        error('plotDipoles requires at least two inputs: xyzMatrix and rgbMatrix.');
    end
    if nargin < 3
        mriData = []; % Default to no MRI if not provided
    elseif ~isempty(mriData) && (~isstruct(mriData) && ~ischar(mriData))
        warning('mriData is not a valid structure or file path. Proceeding without MRI.');
        mriData = [];
    end

    % Build the dipole data structure and RGB color cells
    dipModel = struct('posxyz', [], 'momxyz', [], 'rv', []);
    rgbCells = cell(size(xyzMatrix, 1), 1);
    for dipIdx = 1:size(xyzMatrix, 1)
        dipModel(dipIdx).posxyz = xyzMatrix(dipIdx, :);
        dipModel(dipIdx).momxyz = [0 0 0];
        dipModel(dipIdx).rv = 0;
        rgbCells{dipIdx} = rgbMatrix(dipIdx, :);
    end

    % Plot dipoles (spheres only) into the current figure/axes
    % If no figure exists, the caller should create one
    try
        if isempty(mriData)
            dipPlotHandle = dipplot(dipModel, 'coordformat', 'MNI', 'color', rgbCells, ...
                                   'gui', 'off', 'dipolesize', 30, 'view', [45 35], ...
                                   'spheres', 'on', 'projlines', 'off');
        else
            dipPlotHandle = dipplot(dipModel, 'coordformat', 'MNI', 'color', rgbCells, ...
                                   'gui', 'off', 'dipolesize', 30, 'view', [45 35], ...
                                   'spheres', 'on', 'projlines', 'off', 'mri', mriData);
        end
    catch ME
        warning(ME.identifier, 'Failed to plot dipoles: %s. Check mriData format or dipplot compatibility.', ME.message);
        dipPlotHandle = [];
    end
end