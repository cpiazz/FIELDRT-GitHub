function [planC, jaccard3Doutput, jaccard2Doutput, volumestat_GS, volumestat_US, COM_GS, COM_US, volume_Ratio, OARs_finalmask, OARsdummy_finalmask] = FIELDRTAnalysis(planC, Case, FIELDRTGSCases)
% FIELDRTAnalysis
% Computes the analysis
%
% 02/10/2017
%
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

% import org.dcm4che2.io.*
% import org.dcm4che2.data.*

warning off

%%  Computing inters_contour, union_contour, GS_conotur, user_contour, some conformity indexes
%% and the possible intersection between the US structures and the OARs
% inters_contour = intersection between GS and user contour;
% union_contour = union between GS and user contour;
% GS_conotur = difference between GS contour and intesected contour;
% user_contour = difference between user contour and intesected contour;
% conformity indexes = volume, COM, Volume Ratio, 3D Jaccard, 2D Jaccard;

indexS = planC{end};
numStruct = indexS.structures;
numofStruct = size(planC{1, numStruct}, 2);
scanNum = 1;
OARs_finalmask = []; % Mask for the OARs
OARsdummy_finalmask = []; % Mask for the dummy OARs
% To be changed based on the case analysed
% Case = 1;

% if Case == 11 % Oesophagus Case1
%     struct_name_case = {'GTV_GS', 'CTVA_GS', 'CTVB_GS', 'CTVC_GS', 'PTV_GS'};    
% else
%     
% end
% 
% if Case == 31 % Prostate Case1
%     struct_name_case = {'CTVp_GS', 'CTVpsv_GS'};    
% else
%     
% end

struct_name_case = FIELDRTGSCases(Case(1,1)).Structures(Case(1,2)).Cases;

waitingbar_index = 0;

% Waiting bar
wb = waitbar(waitingbar_index, 'Computing the analyses ...');


% Computing new contours for each structure
for k = 1 : size(struct_name_case, 2) % (numofStruct / 2)
    
    % Union and intersection for all structures in planC
    struct_name_temp = FIELDRTGSCases(Case(1,1)).StructOverunder(Case(1,2)).Cases(k); % planC{1, 4}(k).structureName;  % Name of the structure to analyse. You need this for converting back the back to RasterSegments
    
    %     structRefNum = planC{1, 4}(k).structureName;  % Name of the GS structure to analyse
    %     structEvalNum = planC{1, 4}(k + (numofStruct / 2)).structureName;  % Name of the user structure to analyse
    %   Getting the number of the structures segmented by the user to analyse
    structRefNum  = getStructNum([struct_name_temp{1, 1} '_GS'], planC, indexS);
    structEvalNum = getStructNum(struct_name_temp, planC, indexS);
   
    % Union between structures
    % [rasterSegs_union, unionVolume] = FIELDRTcalcUnionVol(structEvalNum, structRefNum, planC, planC, struct_name_temp);
    [rasterSegs_union, unionMask] = FIELDRTcalcUnionVol(planC{1, 4}(structRefNum).rasterSegments, planC{1, 4}(structEvalNum).rasterSegments, scanNum, planC);  % And if there is more than one scan?
    
    % Copying the information for the new structure
    planC{1, 4}(size(planC{1, numStruct}, 2) + 1) = planC{1, 4}(structRefNum);
    numstructtoupdate = size(planC{1, numStruct}, 2);
    
    % Adding the new modified unified structure
    planC{1, 4}(size(planC{1, numStruct}, 2)).structureName = [struct_name_temp{1, 1} '_union']; % Name of the new modified structure
    planC{1, 4}(size(planC{1, numStruct}, 2)).rasterSegments = rasterSegs_union; % Adding the name of the rasterSegments
    planC{1, 4}(size(planC{1, numStruct}, 2)).contour = (FIELDRTrasterToPoly(rasterSegs_union, 1, planC))'; % Adding the contours
    planC = updateStructureMatrices(planC, size(planC{1, numStruct}, 2)); %% Updating uniformized data for changes in structure contours
    
    % Intersection between structures
    % [rasterSegs_inters, intersectVolume] = FIELDRTcalcIntersectionVol(structEvalNum, structRefNum, planC, planC, struct_name_temp);
    [rasterSegs_inters, intersectMask] = FIELDRTcalcIntersectionVol(planC{1, 4}(structRefNum).rasterSegments, planC{1, 4}(structEvalNum).rasterSegments, scanNum, planC);  % And if there is more than one scan?
    
    % Copying the information for the new structure
    planC{1, 4}(size(planC{1, numStruct}, 2) + 1) = planC{1, 4}(structRefNum);
    numstructtoupdate = [numstructtoupdate; size(planC{1, numStruct}, 2)];
    
    % Updating the waiting bar
    waitbar(((waitingbar_index + 1)/ (size(struct_name_case, 2) * 5)), wb);
    
    % Adding the new modified intersected structure
    planC{1, 4}(size(planC{1, numStruct}, 2)).structureName = [struct_name_temp{1, 1} '_inters']; % Name of the new modified structure
    planC{1, 4}(size(planC{1, numStruct}, 2)).rasterSegments = rasterSegs_inters; % Adding the name of the rasterSegments
    planC{1, 4}(size(planC{1, numStruct}, 2)).contour = (FIELDRTrasterToPoly(rasterSegs_inters, 1, planC))'; % Adding the contours
    planC = updateStructureMatrices(planC, size(planC{1, numStruct}, 2));  %% Updating uniformized data for changes in structure contours
    
    % Difference between GS contour and intesected contour
    % [rasterSegs_diff_GS] = FIELDRTcalcDifferenceVol(size(planC{1, numStruct}, 2), structRefNum, planC, planC, struct_name_temp);
    [rasterSegs_diff_GS, diffMask_GS] = FIELDRTcalcDifferenceVol(planC{1, 4}(structRefNum).rasterSegments, planC{1, 4}(structEvalNum).rasterSegments, scanNum, planC);
    
    % Copying the information for the new structure
    planC{1, 4}(size(planC{1, numStruct}, 2) + 1) = planC{1, 4}(structRefNum);
    numstructtoupdate = [numstructtoupdate; size(planC{1, numStruct}, 2)];
    
    % Updating the waiting bar
    waitbar(((waitingbar_index + 2)/ (size(struct_name_case, 2) * 5)), wb);
    
    % Adding the new modified intersected structure
    planC{1, 4}(size(planC{1, numStruct}, 2)).structureName = [struct_name_temp{1, 1} '_GSdiffinters']; % Name of the new modified structure
    planC{1, 4}(size(planC{1, numStruct}, 2)).rasterSegments = rasterSegs_diff_GS; % Adding the name of the rasterSegments
    planC{1, 4}(size(planC{1, numStruct}, 2)).contour = (FIELDRTrasterToPoly(rasterSegs_diff_GS, 1, planC))'; % Adding the contours
    planC = updateStructureMatrices(planC, size(planC{1, numStruct}, 2));  %% Updating uniformized data for changes in structure contours
    
    % Difference between user contour and intesected contour
    % [rasterSegs_inters_US] = FIELDRTcalcDifferenceVol(size(planC{1, numStruct}, 2) - 1, structEvalNum, planC, planC, struct_name_temp);
    [rasterSegs_diff_US, diffMask_US] = FIELDRTcalcDifferenceVol(planC{1, 4}(structEvalNum).rasterSegments, planC{1, 4}(structRefNum).rasterSegments, scanNum, planC);
    
    % Copying the information for the new structure
    planC{1, 4}(size(planC{1, numStruct}, 2) + 1) = planC{1, 4}(structEvalNum);
    numstructtoupdate = [numstructtoupdate; size(planC{1, numStruct}, 2)];
    
    % Updating the waiting bar
    waitbar(((waitingbar_index + 3)/ (size(struct_name_case, 2) * 5)), wb);
    
    % Adding the new modified intersected structure
    planC{1, 4}(size(planC{1, numStruct}, 2)).structureName = [struct_name_temp{1, 1} '_USdiffinters']; % Name of the new modified structure
    planC{1, 4}(size(planC{1, numStruct}, 2)).rasterSegments = rasterSegs_diff_US; % Adding the name of the rasterSegments
    planC{1, 4}(size(planC{1, numStruct}, 2)).contour = (FIELDRTrasterToPoly(rasterSegs_diff_US, 1, planC))'; % Adding the contours
    planC = updateStructureMatrices(planC, size(planC{1, numStruct}, 2));  %% Updating uniformized data for changes in structure contours
    
    %     % Updating the matrices
    %     for index3 = 1 : size(numstructtoupdate, 1)
    %         planC = updateStructureMatrices(planC, numstructtoupdate(index3));  %% Updating uniformized data for changes in structure contours
    %     end
    
    % Updating the waiting bar
    waitbar(((waitingbar_index + 4)/ (size(struct_name_case, 2) * 5)), wb);
    
    %%  Computing conformity indexes (Volume, COM, Volume Ratio, 3D Jaccard, 2D Jaccard, Mean Distance to Conformity (MDC))
    volumestat_GS{k, 1} = FIELDRTVolumecomputation(structRefNum, planC); % GS
    volumestat_US{k, 1} = FIELDRTVolumecomputation(structEvalNum, planC); % US
    
    volume_Ratio{k, 1} = volumestat_GS{k, 1}.vol / volumestat_US{k, 1}.vol;
    
    COM_GS{k, 1} = FIELDRTCOMcomputation(structRefNum, planC); % GS
    COM_US{k, 1} = FIELDRTCOMcomputation(structEvalNum, planC); % US
    
    jaccard3Doutput{k, 1} = nnz(intersectMask) / nnz(unionMask);
    
    for indexjac = 1 : size(intersectMask, 3)
        jaccard2Doutput{k, indexjac} = nnz(intersectMask(:, :, indexjac)) / nnz(unionMask(:, :, indexjac));
    end    
    
    %% Computing the possible intersection between US structures and the shrink OARs.
    
%     if Case == 11 % Computed only for GTV, CTVB and PTV(removed)
%         if (strcmp(struct_name_case{k}, 'GTV_GS') || strcmp(struct_name_case{k}, 'CTVB_GS')) % || strcmp(struct_name_case{k}, 'PTV_GS'))
%             struct_name_temp = struct_name_case{k};
%             structEvalNum = getStructNum(struct_name_temp(1 : end - 3), planC, indexS);
%             
%             OARs_struct_name_case = {'Vertebra'; 'Aorta'; 'Right lung'; 'Pericardium/great vessels'; 'Liver'; 'Stomach'; 'Azygous vein'; 'Left main bronchus'; 'Left lung'};
%             
%             sizeOARstemp = size(OARsdummy_finalmask, 1) + 1;
%             for indexOARs = 1 : size(OARs_struct_name_case, 1)
%                 % Intersection between the US outline and the dummyOARs
%                 %   Getting the number of the structures segmented by the user to analyse
%                 structRefNum  = getStructNum([OARs_struct_name_case{indexOARs} '_shrink_' struct_name_case{k}], planC, indexS);
%                 structEvalNum = getStructNum(struct_name_temp(1 : end - 3), planC, indexS);
%                 
%                 % intersection between the US contour and the dummy OARs
%                 [rasterSegs_inters, intersectMask_dummyOARs] = FIELDRTcalcIntersectionVol(planC{1, 4}(structRefNum).rasterSegments, planC{1, 4}(structEvalNum).rasterSegments, scanNum, planC);  % And if there is more than one scan?
%                 OARsdummy_finalmask{sizeOARstemp, indexOARs} = intersectMask_dummyOARs;
%                 
%                 % Intersection between the US outline and the OARs (not
%                 % needed)
% %                 structRefNum  = getStructNum([OARs_struct_name_case{indexOARs} '_GS'], planC, indexS);
% %                 
% %                 % intersection between the US contour and the OARs (used in the
% %                 % viewer and the pdf to obtain the centroid
% %                 [rasterSegs_inters, intersectMask_OARs] = FIELDRTcalcIntersectionVol(planC{1, 4}(structRefNum).rasterSegments, planC{1, 4}(structEvalNum).rasterSegments, scanNum, planC);  % And if there is more than one scan?
% %                 OARs_finalmask{sizeOARstemp, indexOARs} = intersectMask_OARs;
%             end
%         end
%         
%     end
%     
%     if Case == 31 % Prostate Case1
%        
%             struct_name_temp = struct_name_case{k};
%             structEvalNum = getStructNum(struct_name_temp(1 : end - 3), planC, indexS);
%             
%             OARs_struct_name_case = {'Bladder'; 'Bowel'; 'Lt FemHead'; 'Penile bulb'; 'Rectum'; 'Rt FemHead'};
%             
%             
%             sizeOARstemp = size(OARsdummy_finalmask, 1) + 1;
%             for indexOARs = 1 : size(OARs_struct_name_case, 1)
%                 % Intersection between the US outline and the dummyOARs
%                 %   Getting the number of the structures segmented by the user to analyse
%                 structRefNum  = getStructNum([OARs_struct_name_case{indexOARs} '_shrink_' struct_name_case{k}], planC, indexS);
%                 structEvalNum = getStructNum(struct_name_temp(1 : end - 3), planC, indexS);
%                 
%                 % intersection between the US contour and the dummy OARs
%                 [rasterSegs_inters, intersectMask_dummyOARs] = FIELDRTcalcIntersectionVol(planC{1, 4}(structRefNum).rasterSegments, planC{1, 4}(structEvalNum).rasterSegments, scanNum, planC);  % And if there is more than one scan?
%                 OARsdummy_finalmask{sizeOARstemp, indexOARs} = intersectMask_dummyOARs;
%                 
%                 % Intersection between the US outline and the OARs (not
%                 % needed)
%                 %                 structRefNum  = getStructNum([OARs_struct_name_case{indexOARs} '_GS'], planC, indexS);
%                 %
%                 %                 % intersection between the US contour and the OARs (used in the
%                 %                 % viewer and the pdf to obtain the centroid
%                 %                 [rasterSegs_inters, intersectMask_OARs] = FIELDRTcalcIntersectionVol(planC{1, 4}(structRefNum).rasterSegments, planC{1, 4}(structEvalNum).rasterSegments, scanNum, planC);  % And if there is more than one scan?
%                 %                 OARs_finalmask{sizeOARstemp, indexOARs} = intersectMask_OARs;
%             end        
%     end
    
    struct_name_temp_OARs = FIELDRTGSCases(Case(1,1)).StructOARs(Case(1,2)).Cases;
    OARs_struct_name_case = FIELDRTGSCases(Case(1,1)).OARs(Case(1,2)).Cases;
    
    % To check if the current US volume has to be compared with OARs and if it is so the comparison is performed    
    for indexOARs = 1 : size(struct_name_temp_OARs, 2)
        if strcmp(struct_name_temp{1, 1}, struct_name_temp_OARs{1, indexOARs})       
            
            sizeOARstemp = size(OARsdummy_finalmask, 1) + 1;
            for indexOARs = 1 : size(OARs_struct_name_case, 2)
                % Intersection between the US outline and the dummyOARs
                %   Getting the number of the structures segmented by the user to analyse
                structRefNum  = getStructNum([OARs_struct_name_case{indexOARs} '_shrink_' struct_name_temp{1, 1} '_GS'], planC, indexS);
                structEvalNum = getStructNum(struct_name_temp{1, 1}, planC, indexS);
                
                % intersection between the US contour and the dummy OARs
                [rasterSegs_inters, intersectMask_dummyOARs] = FIELDRTcalcIntersectionVol(planC{1, 4}(structRefNum).rasterSegments, planC{1, 4}(structEvalNum).rasterSegments, scanNum, planC);  % And if there is more than one scan?
                OARsdummy_finalmask{sizeOARstemp, indexOARs} = intersectMask_dummyOARs;
                
                % Intersection between the US outline and the OARs (not
                % needed)
                %                 structRefNum  = getStructNum([OARs_struct_name_case{indexOARs} '_GS'], planC, indexS);
                %
                %                 % intersection between the US contour and the OARs (used in the
                %                 % viewer and the pdf to obtain the centroid
                %                 [rasterSegs_inters, intersectMask_OARs] = FIELDRTcalcIntersectionVol(planC{1, 4}(structRefNum).rasterSegments, planC{1, 4}(structEvalNum).rasterSegments, scanNum, planC);  % And if there is more than one scan?
                %                 OARs_finalmask{sizeOARstemp, indexOARs} = intersectMask_OARs;
            end
        end
    end 
    
    
    % Updating the waiting bar
    waitbar(((waitingbar_index + 5)/ (size(struct_name_case, 2) * 5)), wb);
    waitingbar_index = waitingbar_index + 5;
end

% Close waiting bar
close(wb);

end
