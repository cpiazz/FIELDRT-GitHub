function [sag_masks, cor_masks] = FIELDRTconvertstructsview(planC)
% Convert the axial contours in sagittal and coronal

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
% Concetta Piazzese March 2018

indexS = planC{end};
numScan = indexS.scan;
numStruct = indexS.structures;
numofStruct = size(planC{1, numStruct}, 2);
scanNum = 1; % To change if there are more scans!!!!
masks = [];

% Computing a mask for each structure in the sagittal and coronal view

% First a a mask with the axial contour will be obtained. 
% The mask will be then rotated and saved
%

%% The sagittal and coronal mask of the GS should be created when setting up the GS!!!!

for k = 1 : numofStruct % (numofStruct / 2)
    
    struct_name_temp = planC{1, numStruct}(k).structureName; % planC{1, 4}(k).structureName;  % Name of the structure to analyse. You need this for converting back the back to RasterSegments
    
    %   Getting the number of the structures to copy
    structRefNum  = getStructNum(struct_name_temp, planC, indexS);
    
    rasterSegs_struct = [];
    
    if isempty(planC{1, numStruct}(structRefNum).rasterSegments) % if there are no raster segments
        mask_struct = zeros(size(planC{numScan}(scanNum).scanArray));
        
    else % if there are raster segments
        rasterSegs_struct = planC{1, numStruct}(structRefNum).rasterSegments;
        
        scanNum = getStructureAssociatedScan(structRefNum, planC);
        
        %sort input rasterSegments by CTSliceValue
        rasterSegs_struct = sortrows(rasterSegs_struct, 6);
        
        rasterIndices = [];
        mask_struct = [];
        for i = 1 : size(planC{numScan}(scanNum).scanArray, 3) % length(slicesToCalculate)
            sliceNum = i; % slicesToCalculate(i);
            rasterIndices = find(rasterSegs_struct(:,6) == sliceNum);
            mask_struct(:, :, i) = rasterToMask(rasterSegs_struct(rasterIndices,:), scanNum, planC);
        end        
    end
    % Sagittal
    mask_s = [];
    % sag_mask_bound = [];
    mask_s = squeeze(permute(mask_struct, [3 1 2 4])); % Sagittal view image
    % sag_mask_bound = bwboundaries(sag_mask(:, :, handles.currentSlice_sag));
    
    sag_masks(k).structureName = [struct_name_temp(1 : end - 2) 'sag']; % Name of the new modified structure
    sag_masks(k).mask = mask_s; % Adding the mask
    
    %     Xrange = handles.sagittalxlim(2) - handles.sagittalxlim(1);
    %     Yrange = handles.sagittalylim(2) - handles.sagittalylim(1);
    
    % Coronal
    mask_c = [];
    % cor_mask_bound = [];
    mask_c = squeeze(permute(mask_struct, [3 2 1 4])); % Coronal view image
    %     cor_mask_bound = bwboundaries(cor_mask(:, :, handles.currentSlice_sag));
    %
    cor_masks(k).structureName = [struct_name_temp(1 : end - 2) 'cor']; % Name of the new modified structure
    cor_masks(k).mask = mask_c; % Adding the mask
    %
    %     Xrange = handles.sagittalxlim(2) - handles.sagittalxlim(1);
    %     Yrange = handles.sagittalylim(2) - handles.sagittalylim(1);
    
    
end





