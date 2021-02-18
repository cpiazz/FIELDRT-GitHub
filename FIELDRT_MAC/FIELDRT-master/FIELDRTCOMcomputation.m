function [output] = FIELDRTCOMcomputation(structNum, planC)
% Calculates COM of a structure

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
% Concetta Piazzese October 2018

indexS = planC{end};

indexS = planC{end};

%Get scan associated with structure.
scanSet = getStructureAssociatedScan(structNum, planC);

%Get mask of structure.
uniformStr = getUniformStr(structNum, planC);

%Get r,c,s of voxels inside uniformStr.
[r,c,s] = find3d(uniformStr);
nVoxInStruct = length(r);
clear uniformStr;

%Get scan's original x,y,z coordinates, and it's transformation matrix.
if isempty(planC{indexS.scan}(scanSet).uniformScanInfo)
    planC = setUniformizedData(planC);
end
[xV, yV, zV] = getUniformScanXYZVals(planC{indexS.scan}(scanSet));
transM = getTransM(planC{indexS.scan}(scanSet), planC);
voxVol = abs(xV(2)-xV(1)) * abs(yV(2)-yV(1)) * abs(zV(2)-zV(1));

%Get the x,y,z coords of points in the structure, and transform.
structXV = xV(c); clear c xV;
structYV = yV(r); clear r yV;
structZV = zV(s); clear s zV;
[xT, yT, zT] = applyTransM(transM, structXV, structYV, structZV);
clear structXV structYV structZV

%Find the bounding box around the structure.
minX = min(xT);
minY = min(yT);
minZ = min(zT);
maxX = max(xT);
maxY = max(yT);
maxZ = max(zT);

yCOM = mean(yT);
xCOM = mean(xT);
zCOM = mean(zT);

output.COM = [xCOM, yCOM, zCOM];

end