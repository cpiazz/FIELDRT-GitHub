function contour = FIELDRTrasterToPoly(rasterSegments, scanNum, planC)
% Modified version of rasterToPoly funcion of CERR. The difference is that
% if the rasterSegments matrix is empty the contour matrix of planC is not
% left empty but is populated with the right (empty) field.

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
% Concetta Piazzese March 2020

if isempty(rasterSegments)
    indexS = planC{end};
    [xSize,ySize,zSize] = size(getScanArray(planC{indexS.scan}(scanNum)));
    
    for sliceNum = 1:zSize
        contour(sliceNum).segments.points = [];
    end
    if length(contour) ~= zSize
        contour(zSize).segments.points = [];
    end
else
    
    indexS = planC{end};
    [xSize,ySize,zSize] = size(getScanArray(planC{indexS.scan}(scanNum)));
    
    %sort input rasterSegments by CTSliceValue
    rasterSegments = sortrows(rasterSegments, 6);
    
    %get list of CTSlices to iterate over.
    slicesToCalculate = unique(rasterSegments(:,6));
    
    for sliceNum = 1:zSize
        if any(slicesToCalculate==sliceNum)
            rasterIndices = find(rasterSegments(:,6) == sliceNum);
            mask = rasterToMask(rasterSegments(rasterIndices,:), scanNum, planC);
            contour(sliceNum) = maskToPoly(mask, sliceNum, scanNum, planC);
        else
            contour(sliceNum).segments.points = [];
        end
    end
    if length(contour) ~= zSize
        contour(zSize).segments.points = [];
    end
end