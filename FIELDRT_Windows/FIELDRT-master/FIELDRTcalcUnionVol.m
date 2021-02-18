function [rasterSegs_union, unionMask] = FIELDRTcalcUnionVol(rasterSegs1, rasterSegs2, scanNum, planC)
% function [unionVolume, intersectVolume, unionStr3M, intersectStr3M] = calcInterUnionIntersectionVol_CP(structEvalNum,structRefNum,planCeval,planCref)
% Modified version of structUnion function contained in CERR

% Calculates union volume between two structures from
% different CERR containers.
% Assumes that structures are associated with the same scan.

% Copyright 2010, Joseph O. Deasy, on behalf of the CERR development team.
% Modified by Concetta Piazzese March 2018

if isempty(rasterSegs1) || isempty(rasterSegs2) 
    rasterSegs_union = [];
    unionMask = logical(zeros(size(planC{1, 3}.scanArray)));
else
    indexS = planC{end};
    rasterSegs_union = [];
    
    % %get list of CTSlices to iterate over.
    % slices1 = unique(rasterSegs1(:,6));
    % slices2 = unique(rasterSegs2(:,6));
    %
    % %for union, need to worry about slices where either structure has segments.
    % slicesToCalculate = union(slices1, slices2);
    
    
    % sort input rasterSegments by CTSliceValue
    rasterSegs1 = sortrows(rasterSegs1, 6);
    rasterSegs2 = sortrows(rasterSegs2, 6);
    
    %for each slice we are calculating on, create a mask for each structure and
    %intersect them. Then convert from that mask to raster segments.
    for i = 1 : size(planC{1, 3}.scanArray, 3) % length(slicesToCalculate)
        sliceNum = i; % slicesToCalculate(i);
        rasterIndices = find(rasterSegs1(:,6) == sliceNum);
        mask1 = rasterToMask(rasterSegs1(rasterIndices,:), scanNum, planC);
        rasterIndices = find(rasterSegs2(:,6) == sliceNum);
        mask2 = rasterToMask(rasterSegs2(rasterIndices,:), scanNum, planC);
        unionMasktemp = mask1 | mask2;
        rasterSegs_union = [rasterSegs_union;maskToRaster(unionMasktemp, sliceNum, scanNum, planC)];
        unionMask(:, :, i) = unionMasktemp;
    end
end
end