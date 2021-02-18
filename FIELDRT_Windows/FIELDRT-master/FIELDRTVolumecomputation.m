function [output] = FIELDRTVolumecomputation(structNum, planC)
% Calculates volume of a structure

% FIELDRT is distributed under the terms of the Lesser GNU Public License.
%
%     This version of FIELDRT is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
% FIELDRT is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
% without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
% See the GNU General Public License for more details.
%
% Author
% Concetta Piazzese September 2018

indexS = planC{end};

% Get mask of structure.
mask_struct = double(getUniformStr(structNum, planC));

% Get scan associated with structure.
scanNum    = getStructureAssociatedScan(structNum, planC);

% Record voxel spacing [in mm] (needed for morphology features)
[~, spacing, ~] = getScanOriginSpacing(planC{indexS.scan}(scanNum));

voxelDim = spacing*10; 

% Computing the volume
volumeResult = (sum(mask_struct(:)).*prod(voxelDim));

output.vol = volumeResult/1000;    

end