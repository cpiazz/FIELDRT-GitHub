function FIELDRTplotContours(handles, index_selected)
% Plot the contours

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

% Retrieving some information
if ~isempty(handles.planC{handles.scanfieldnum}.scanInfo(1).xOffset)
    xCTOffset = handles.planC{handles.scanfieldnum}.scanInfo(1).xOffset;
    yCTOffset = handles.planC{handles.scanfieldnum}.scanInfo(1).yOffset;
else
    xCTOffset = 0;
    yCTOffset = 0;
end


%xAAPMShifted -- x coordinates in AAPM system assuming midpoint may be shifted from center
%yAAPMShifted -- y coordinates in AAPM system assuming midpoint may be shifted from center
%xOffset -- x offset from center point
%yOffset -- y offset from center point
%ImageWidth -- a single number (implying a square image) or a length-2 vector,
%              giving the numbers of rows first then columns.
%voxelSizeV -- length-2 vector giving [length on the y-side, length on the x-side]

scale =  handles.planC{handles.scanfieldnum}.scanInfo(handles.currentSlice_ax).grid1Units;

imageWidth(1, 1) =  handles.planC{handles.scanfieldnum}.scanInfo(handles.currentSlice_ax).sizeOfDimension1;
imageWidth(2, 1) =  handles.planC{handles.scanfieldnum}.scanInfo(handles.currentSlice_ax).sizeOfDimension2;

% Deleting children objects corresponding to the contours. Doing it this
% way because set(handles.axes3.Children() command is not working properly
% and the green contour is not the last one appearing.
childsize = size(handles.axes3.Children, 1);
if childsize > 1
    for indexchild = 1 : childsize - 1
        delete(handles.axes3.Children(1));
    end
end

childsize = size(handles.axes5.Children, 1);
if childsize > 1
    for indexchild = 1 : childsize - 1
        delete(handles.axes5.Children(1));
    end
end

childsize = size(handles.axes6.Children, 1);
if childsize > 1
    for indexchild = 1 : childsize - 1
        delete(handles.axes6.Children(1));
    end
end

indexS = handles.planC{end};

% Deleting all annotations in the image
% delete(findall(gcf,'type','annotation'))

% if isempty(index_selected)
%     
% else
%     % To be changed based on the case analysed
%     % Case = 1;
%     indexS = handles.planC{end};
%     
%     if handles.Case == 11
%         struct_name_case = {'GTV', 'CTVA', 'CTVB', 'CTVC', 'PTV'};
%         struct_name_temp = struct_name_case{index_selected};
%         % Retrieving information about the selected structure
%         GSsstructure_original  = getStructNum([struct_name_temp '_GS'], handles.planC, indexS);
%         USsstructure_original  = getStructNum([struct_name_temp], handles.planC, indexS);
%         GSsstructure  = getStructNum([struct_name_temp '_GSdiffinters'], handles.planC, indexS);
%         USsstructure  = getStructNum([struct_name_temp '_USdiffinters'], handles.planC, indexS);
%         Intersectedstructures  = getStructNum([struct_name_temp '_inters'], handles.planC, indexS);
%         GSsstructuremax  = getStructNum([struct_name_temp '_max_GS'], handles.planC, indexS);
%         GSsstructuremin  = getStructNum([struct_name_temp '_min_GS'], handles.planC, indexS);
%     end
%     
%     if handles.Case == 31
%         struct_name_case = {'CTVp', 'CTVpsv'};
%         struct_name_temp = struct_name_case{index_selected};
%         % Retrieving information about the selected structure
%         GSsstructure_original  = getStructNum([struct_name_temp '_GS'], handles.planC, indexS);
%         USsstructure_original  = getStructNum([struct_name_temp], handles.planC, indexS);
%         GSsstructure  = getStructNum([struct_name_temp '_GSdiffinters'], handles.planC, indexS);
%         USsstructure  = getStructNum([struct_name_temp '_USdiffinters'], handles.planC, indexS);
%         Intersectedstructures  = getStructNum([struct_name_temp '_inters'], handles.planC, indexS);
%         GSsstructuremax  = getStructNum([struct_name_temp '_max_GS'], handles.planC, indexS);
%         GSsstructuremin  = getStructNum([struct_name_temp '_min_GS'], handles.planC, indexS);
%     end    
% end



% color coding:
% GS's structure --> blue
% User's structure --> red
% Intersected structures --> green

hhl = [];
labelhhl = [];
hhlmask = []; % For OARs shrink
labelhhlmask = []; % For OARs shrink

% Over/under contoured regions --> Case1
if handles.overunder.Value == 1 && handles.oars.Value == 0 && handles.accreg.Value == 0 && handles.refvolume.Value == 0
    n = 1;
end

% Oars --> Case2
if handles.overunder.Value == 0 && handles.oars.Value == 1 && handles.accreg.Value == 0 && handles.refvolume.Value == 0
    n = 2;
end

% Min/max contoured regions --> Case3
if handles.overunder.Value == 0 && handles.oars.Value == 0 && handles.accreg.Value == 1 && handles.refvolume.Value == 0
    n = 3;
end

% reference volume --> Case4
if handles.overunder.Value == 0 && handles.oars.Value == 0 && handles.accreg.Value == 0 && handles.refvolume.Value == 1
    n = 4;
end

% every toogle button deselected --> Case5
if handles.overunder.Value == 0 && handles.oars.Value == 0 && handles.accreg.Value == 0 && handles.refvolume.Value == 0
    n = 5;
end


%% Axial view

% Setting the axes where visualize the contours
axes(handles.axes3)

% Information about the switch
% case1 = over/under contoured regions
% case2 = OARs
% case3 = max/min
% case4 = GS/user contour (no intersection)

% Emptying the vector
struct_name_case = [];

switch n
    case 1
        if isempty(index_selected) % When no structure is selected
            % legend(handles.axes3,'hide')
            
        else
            % Retrieving information about the structures to show
            struct_name_case = handles.FIELDRTGSCases(handles.Case(1,1)).StructOverunder(handles.Case(1,2)).Cases;
            GSsstructure_original  = getStructNum([struct_name_case{1, index_selected} '_GS'], handles.planC, indexS);
            USsstructure_original  = getStructNum(struct_name_case{1, index_selected}, handles.planC, indexS);
            GSsstructure  = getStructNum([struct_name_case{1, index_selected} '_GSdiffinters'], handles.planC, indexS);
            USsstructure  = getStructNum([struct_name_case{1, index_selected} '_USdiffinters'], handles.planC, indexS);
            Intersectedstructures  = getStructNum([struct_name_case{1, index_selected} '_inters'], handles.planC, indexS);
            
            %% Plotting Under region contour
            if get(handles.refvolume_checkbox, 'Value') == 1 % plot modified GS, modified US and intersection structure if reference checkbox is activated
                % GS's structure
                if isempty(handles.planC{handles.structfieldnum}(GSsstructure).contour)
                    xCoords = [];
                    yCoords = [];
                    %         set(handles.axes3.Children(1), 'XData', xCoords, 'YData', yCoords, 'Color', [0 1 1], 'linewidth', 1,...
                    %             'visible', 'off')
                    hhl = hhl;
                else
                    hL1 = [];
                    numSegs = length(handles.planC{handles.structfieldnum}(GSsstructure).contour(handles.currentSlice_ax).segments);
                    
                    %     thisax = handles.axes3;
                    %     plot2handle = get(thisax, 'UserData');
                    
                    yCoordsGS = [];
                    xCoordsSG = [];
                    for indexseg = 1 : numSegs
                        pointsM = handles.planC{handles.structfieldnum}(GSsstructure).contour(handles.currentSlice_ax).segments(indexseg).points;
                        if isempty(pointsM) % Where there are no contours to plot
                            xCoords = [];
                            yCoords = [];
                            %                 set(handles.axes3.Children(1), 'XData', [], 'YData', [], 'Color', [0 1 1], 'linewidth', 1,...
                            %                     'visible', 'off')
                            hhl = hhl;
                            
                        else
                            xCoords = pointsM(:,1);
                            yCoords = pointsM(:,2);
                            
                            [rowV, colV] = aapmtom(xCoords/scale, yCoords/scale, xCTOffset/scale, yCTOffset/scale, imageWidth);
                            
                            rowV = [rowV; rowV(1)];
                            colV = [colV; colV(1)];
                            
                            %Update contour plot
                            hold on
                            % hPlot = plot([colV; colV(1)], [rowV; rowV(1)], 'rs-', 'Color', [0 1 1], 'LineStyle', '-');
                            %                 set(handles.axes3.Children(1), 'XData', colV, 'YData', rowV, 'Color', [0 1 1], 'linewidth', 1,...
                            %                     'visible', 'on')
                            hL1 = plot(colV, rowV, 'Color', [1 1 0], 'linewidth', round(get(handles.Slider_zoom, 'Value')));
                            refreshdata(handles.axes3) %, 'Position', [0.1, 0.012])
                            
                            %                     imshow(blue);
                            %                 blue = cat(3, zeros(handles.volsize(1 : 2)), zeros(handles.volsize(1 : 2)), ones(handles.volsize(1 : 2)));
                            %                 set(handles.hShow.Parent.Children(7), 'AlphaData', blue);
                            %                 https://blogs.mathworks.com/steve/2009/02/18/image-overlay-using-transparency/
                            
                            %     % GS's structure
                            %     pointsGS = handles.planC{handles.structfieldnum}(4*(index_selected - 1) + 13).contour(handles.currentSlice_ax).segments(indexseg).points;
                            %     yCoordsGS = [yCoordsGS; pointsGS(:,2); pointsGS(1,2)]; % + yT; Rotation component (??)  /Users/concetta/OneDrive - Cardiff
                            %     % University/Concetta_work_cardiff/Projects/Other_projects/CERR-master/CERR_core/Viewers
                            %     % showStructures line 227
                            %     xCoordsSG = [xCoordsSG; pointsGS(:,1); pointsGS(1,1)]; %  + xT; Rotation component (??)  /Users/concetta/OneDrive - Cardiff
                            %     % University/Concetta_work_cardiff/Projects/Other_projects/CERR-master/CERR_core/Viewers
                            %     % showStructures line 227
                        end
                    end
                    
                    % Legend for the contour
                    if isempty(hL1)
                        
                    else
                        hhl = [hhl hL1];
                        labelhhl = [labelhhl {'Under region'}];
                    end
                    
                end
                
                %% Plotting over region contour derived from GS
                if isempty(handles.planC{handles.structfieldnum}(USsstructure).contour)
                    xCoords = [];
                    yCoords = [];
                    %         set(handles.axes3.Children(2), 'XData', xCoords, 'YData', yCoords, 'Color', [1 0 0], 'linewidth', 1,...
                    %             'visible', 'off')
                    hhl = hhl;
                else
                    hL2 = [];
                    numSegs = length(handles.planC{handles.structfieldnum}(USsstructure).contour(handles.currentSlice_ax).segments);
                    
                    yCoordsGS = [];
                    xCoordsSG = [];
                    rowV = [];
                    colV = [];
                    for indexseg = 1 : numSegs
                        pointsM = handles.planC{handles.structfieldnum}(USsstructure).contour(handles.currentSlice_ax).segments(indexseg).points;
                        
                        if isempty(pointsM) % Where there are no contours to plot
                            xCoords = [];
                            yCoords = [];
                            %                 set(handles.axes3.Children(2), 'XData', xCoords, 'YData', yCoords, 'Color', [1 0 0], 'linewidth', 1,...
                            %                     'visible', 'off')
                            hhl = hhl;
                            
                        else
                            xCoords = pointsM(:,1);
                            yCoords = pointsM(:,2);
                            
                            [rowV, colV] = aapmtom(xCoords/scale, yCoords/scale, xCTOffset/scale, yCTOffset/scale, imageWidth);
                            
                            rowV = [rowV; rowV(1)];
                            colV = [colV; colV(1)];
                            
                            %Update contour plot
                            hold on
                            
                            % hPlot = plot([colV; colV(1)], [rowV; rowV(1)], 'rs-', 'Color', [1 0 0], 'LineStyle', '-');
                            %                 set(handles.axes3.Children(2), 'XData', colV, 'YData', rowV, 'Color', [1 0 0], 'linewidth', 1,...
                            %                     'visible', 'on')
                            
                            hL2 = plot(colV, rowV, 'Color', [0 1 1], 'linewidth', round(get(handles.Slider_zoom, 'Value')));
                            refreshdata(handles.axes3)
                            
                            %     % GS's structure
                            %     pointsGS = handles.planC{handles.structfieldnum}(4*(index_selected - 1) + 13).contour(handles.currentSlice_ax).segments(indexseg).points;
                            %     yCoordsGS = [yCoordsGS; pointsGS(:,2); pointsGS(1,2)]; % + yT; Rotation component (??)  /Users/concetta/OneDrive - Cardiff
                            %     % University/Concetta_work_cardiff/Projects/Other_projects/CERR-master/CERR_core/Viewers
                            %     % showStructures line 227
                            %     xCoordsSG = [xCoordsSG; pointsGS(:,1); pointsGS(1,1)]; %  + xT; Rotation component (??)  /Users/concetta/OneDrive - Cardiff
                            %     % University/Concetta_work_cardiff/Projects/Other_projects/CERR-master/CERR_core/Viewers
                            %     % showStructures line 227
                        end
                    end
                    
                    % Legend for the contour
                    if isempty(hL2)
                        
                    else
                        hhl = [hhl hL2];
                        labelhhl = [labelhhl {'Over region'}];
                    end
                end
                
                %% Plotting common region
                % Structures obtained with the intersection of GS and user's structures
                if isempty(handles.planC{handles.structfieldnum}(Intersectedstructures).contour)
                    xCoords = [];
                    yCoords = [];
                    %         set(handles.axes3.Children(3), 'XData', xCoords, 'YData', yCoords, 'Color', [0 1 0], 'linewidth', 1,...
                    %             'visible', 'off')
                    hhl = hhl;
                else
                    hL3 = [];
                    numSegs = length(handles.planC{handles.structfieldnum}(Intersectedstructures).contour(handles.currentSlice_ax).segments);
                    
                    yCoordsGS = [];
                    xCoordsSG = [];
                    rowV = [];
                    colV = [];
                    for indexseg = 1 : numSegs
                        pointsM = handles.planC{handles.structfieldnum}(Intersectedstructures).contour(handles.currentSlice_ax).segments(indexseg).points;
                        
                        if isempty(pointsM) % Where there are no contours to plot
                            xCoords = [];
                            yCoords = [];
                            %                 set(handles.axes3.Children(3), 'XData', xCoords, 'YData', yCoords, 'Color', [0 1 0], 'linewidth', 6,...
                            %                     'visible', 'off')
                            hhl = hhl;
                            
                        else
                            xCoords = pointsM(:,1);
                            yCoords = pointsM(:,2);
                            
                            [rowV, colV] = aapmtom(xCoords/scale, yCoords/scale, xCTOffset/scale, yCTOffset/scale, imageWidth);
                            
                            rowV = [rowV; rowV(1)];
                            colV = [colV; colV(1)];
                            
                            %Update contour plot
                            hold on
                            %             hPlot = plot([colV; colV(1)], [rowV; rowV(1)], 'rs-', 'Color', [0 1 0], 'LineStyle', '-');
                            %                 set(handles.axes3.Children(3), 'XData', colV, 'YData', rowV, 'Color', [0 1 0], 'linewidth', 1,...
                            %                     'visible', 'on')
                            hL3 = plot(colV, rowV, 'Color', [0 1 0], 'linewidth', round(get(handles.Slider_zoom, 'Value')));
                            
                            %     % GS's structure
                            %     pointsGS = handles.planC{handles.structfieldnum}(4*(index_selected - 1) + 13).contour(handles.currentSlice_ax).segments(indexseg).points;
                            %     yCoordsGS = [yCoordsGS; pointsGS(:,2); pointsGS(1,2)]; % + yT; Rotation component (??)  /Users/concetta/OneDrive - Cardiff
                            %     % University/Concetta_work_cardiff/Projects/Other_projects/CERR-master/CERR_core/Viewers
                            %     % showStructures line 227
                            %     xCoordsSG = [xCoordsSG; pointsGS(:,1); pointsGS(1,1)]; %  + xT; Rotation component (??)  /Users/concetta/OneDrive - Cardiff
                            %     % University/Concetta_work_cardiff/Projects/Other_projects/CERR-master/CERR_core/Viewers
                            %     % showStructures line 227
                            refreshdata(handles.axes3)
                            hold on
                        end
                    end
                    
                    % Legend for the contour
                    if isempty(hL3)
                        
                    else
                        hhl = [hhl hL3];
                        labelhhl = [labelhhl {'Common region'}];
                    end
                end
            end
            
            if get(handles.refvolume_checkbox, 'Value') == 0 % plot original user contour if reference checkbox is activated
                %% Plotting over region contour derived from US
                struct_name_temp = struct_name_case{index_selected};
                % Retrieving information about the selected structure
                USsstructure  = getStructNum(struct_name_temp, handles.planC, indexS);
                
                % User's structure
                if isempty(handles.planC{handles.structfieldnum}(USsstructure).contour)
                    xCoords = [];
                    yCoords = [];
                    %         set(handles.axes3.Children(2), 'XData', xCoords, 'YData', yCoords, 'Color', [1 0 0], 'linewidth', 1,...
                    %             'visible', 'off')
                    hhl = hhl;
                else
                    hL2 = [];
                    numSegs = length(handles.planC{handles.structfieldnum}(USsstructure).contour(handles.currentSlice_ax).segments);
                    
                    yCoordsGS = [];
                    xCoordsSG = [];
                    rowV = [];
                    colV = [];
                    for indexseg = 1 : numSegs
                        pointsM = handles.planC{handles.structfieldnum}(USsstructure).contour(handles.currentSlice_ax).segments(indexseg).points;
                        
                        if isempty(pointsM) % Where there are no contours to plot
                            xCoords = [];
                            yCoords = [];
                            %                 set(handles.axes3.Children(2), 'XData', xCoords, 'YData', yCoords, 'Color', [1 0 0], 'linewidth', 1,...
                            %                     'visible', 'off')
                            hhl = hhl;
                            
                        else
                            xCoords = pointsM(:,1);
                            yCoords = pointsM(:,2);
                            
                            [rowV, colV] = aapmtom(xCoords/scale, yCoords/scale, xCTOffset/scale, yCTOffset/scale, imageWidth);
                            
                            rowV = [rowV; rowV(1)];
                            colV = [colV; colV(1)];
                            
                            %Update contour plot
                            hold on
                            
                            % hPlot = plot([colV; colV(1)], [rowV; rowV(1)], 'rs-', 'Color', [1 0 0], 'LineStyle', '-');
                            %                 set(handles.axes3.Children(2), 'XData', colV, 'YData', rowV, 'Color', [1 0 0], 'linewidth', 1,...
                            %                     'visible', 'on')
                            
                            hL2 = plot(colV, rowV, 'Color', [0 1 1], 'linewidth', round(get(handles.Slider_zoom, 'Value')));
                            refreshdata(handles.axes3)
                            
                            %     % GS's structure
                            %     pointsGS = handles.planC{handles.structfieldnum}(4*(index_selected - 1) + 13).contour(handles.currentSlice_ax).segments(indexseg).points;
                            %     yCoordsGS = [yCoordsGS; pointsGS(:,2); pointsGS(1,2)]; % + yT; Rotation component (??)  /Users/concetta/OneDrive - Cardiff
                            %     % University/Concetta_work_cardiff/Projects/Other_projects/CERR-master/CERR_core/Viewers
                            %     % showStructures line 227
                            %     xCoordsSG = [xCoordsSG; pointsGS(:,1); pointsGS(1,1)]; %  + xT; Rotation component (??)  /Users/concetta/OneDrive - Cardiff
                            %     % University/Concetta_work_cardiff/Projects/Other_projects/CERR-master/CERR_core/Viewers
                            %     % showStructures line 227
                        end
                    end
                    
                    % Legend for the contour
                    if isempty(hL2)
                        
                    else
                        hhl = [hhl hL2];
                        labelhhl = [labelhhl {'Over region'}];
                    end
                end
            end
            %             if size(hhl, 2) == 3
            %                 hh = legend(hhl, {'Under region', 'Over region', 'Common region'}, 'FontSize', 12, 'TextColor', 'white', 'Color', 'Black');
            %             end
            %             if size(hhl, 2) == 2
            %                 hh = legend(hhl, {'Under region', 'Over region'}, 'FontSize', 12, 'TextColor', 'white', 'Color', 'Black');
            %             end
            %             if size(hhl, 2) == 1
            %                 if labelhhl(1, 1) == 1
            %                     hh = legend(hhl, {'Under region'}, 'FontSize', 12, 'TextColor', 'white', 'Color', 'Black');
            %                 else
            %                     hh = legend(hhl, {'Over region'}, 'FontSize', 12, 'TextColor', 'white', 'Color', 'Black');
            %                 end
            %             end
            %             if size(hhl, 2) < 1
            %                 legend(handles.axes3,'hide')
            %             end
            if size(hhl, 2) >= 1
                hhllegend = legend(hhl, labelhhl, 'FontSize', 12, 'TextColor', 'white', 'Color', 'Black');
            else
                legend(handles.axes3,'hide')
            end
            
        end
        
    case 2
        if isempty(index_selected) ||  index_selected == 0 % When no structure is selected
            % legend(handles.axes3,'hide')
            % Case = 1; % To be changed based on the analysed case
            indexS = handles.planC{end};
%             % To be changed based on the analysed case
%             if handles.Case == 11
%                 struct_name_case_OARs = {'Vertebra_GS'; 'Aorta_GS'; 'Right lung_GS'; 'Pericardium/great vessels_GS'; ...
%                     'Liver_GS'; 'Stomach_GS'; 'Azygous vein_GS'; 'Left main bronchus_GS'; 'Left lung_GS'};
%                 
%             end
%             
%             if handles.Case == 31
%                 struct_name_case_OARs = {'Bladder_GS'; 'Bowel_GS'; 'Lt FemHead_GS'; 'Penile bulb_GS'; 'Rectum_GS'; 'Rt FemHead_GS'};
%             end
            %% Plotting OARs (even if no structure is selected, OARs have to be showed)
            struct_name_case_OARs = handles.FIELDRTGSCases(handles.Case(1,1)).OARs(handles.Case(1,2)).Cases;
            colorlab = handles.FIELDRTGSCases(handles.Case(1,1)).OARs(handles.Case(1,2)).colorlab;
            
            for indexOARs = 1 : size(struct_name_case_OARs, 2) % Retrieving information about the selected structure
                OARs_struct_name_temp = [struct_name_case_OARs{indexOARs}];
                OARsstructure  = getStructNum([OARs_struct_name_temp '_GS'], handles.planC, indexS);
                
                % Structures obtained with the intersection of GS and user's structures
                if isempty(handles.planC{handles.structfieldnum}(OARsstructure).contour)
                    xCoords = [];
                    yCoords = [];
                    %         set(handles.axes3.Children(3), 'XData', xCoords, 'YData', yCoords, 'Color', [0 1 0], 'linewidth', 1,...
                    %             'visible', 'off')
                    hhl = hhl;
                else
                    hL4 = [];
                    numSegs = length(handles.planC{handles.structfieldnum}(OARsstructure).contour(handles.currentSlice_ax).segments);
                    
                    yCoordsGS = [];
                    xCoordsSG = [];
                    rowV = [];
                    colV = [];
                    for indexseg = 1 : numSegs
                        pointsM = handles.planC{handles.structfieldnum}(OARsstructure).contour(handles.currentSlice_ax).segments(indexseg).points;
                        
                        if isempty(pointsM) % Where there are no contours to plot
                            xCoords = [];
                            yCoords = [];
                            %                 set(handles.axes3.Children(3), 'XData', xCoords, 'YData', yCoords, 'Color', [0 1 0], 'linewidth', 6,...
                            %                     'visible', 'off')
                            hhl = hhl;
                            
                        else
                            xCoords = pointsM(:,1);
                            yCoords = pointsM(:,2);
                            
                            [rowV, colV] = aapmtom(xCoords/scale, yCoords/scale, xCTOffset/scale, yCTOffset/scale, imageWidth);
                            
                            rowV = [rowV; rowV(1)];
                            colV = [colV; colV(1)];
                            
                            %Update contour plot
                            hold on
                            %             hPlot = plot([colV; colV(1)], [rowV; rowV(1)], 'rs-', 'Color', [0 1 0], 'LineStyle', '-');
                            %                 set(handles.axes3.Children(3), 'XData', colV, 'YData', rowV, 'Color', [0 1 0], 'linewidth', 1,...
                            %                     'visible', 'on')
                            
                            
                            
                            hL4 = plot(colV, rowV, 'color', colorlab(indexOARs, :), 'linewidth', round(get(handles.Slider_zoom, 'Value')));
                            
                            %     % GS's structure
                            %     pointsGS = handles.planC{handles.structfieldnum}(4*(index_selected - 1) + 13).contour(handles.currentSlice_ax).segments(indexseg).points;
                            %     yCoordsGS = [yCoordsGS; pointsGS(:,2); pointsGS(1,2)]; % + yT; Rotation component (??)  /Users/concetta/OneDrive - Cardiff
                            %     % University/Concetta_work_cardiff/Projects/Other_projects/CERR-master/CERR_core/Viewers
                            %     % showStructures line 227
                            %     xCoordsSG = [xCoordsSG; pointsGS(:,1); pointsGS(1,1)]; %  + xT; Rotation component (??)  /Users/concetta/OneDrive - Cardiff
                            %     % University/Concetta_work_cardiff/Projects/Other_projects/CERR-master/CERR_core/Viewers
                            %     % showStructures line 227
                            refreshdata(handles.axes3)
                            hold on
                        end
                    end
                    
                    % Legend for the contour
                    % Only for the first segment. We don't need a label for each segment
                    if isempty(hL4)
                        
                    else
                        hhl = [hhl, hL4];
%                         struct_name_temp = struct_name_case_OARs{indexOARs};
%                         labelhhl = [labelhhl {struct_name_temp(1 : end - 3)}];
                        labelhhl = [labelhhl {struct_name_case_OARs{indexOARs}}];
                    end
                end
            end
            
            legend(handles.axes3,'show')
            
            if size(hhl, 2) >= 1
                hhllegend = legend(hhl, labelhhl, 'FontSize', 12, 'TextColor', 'white', 'Color', 'Black');
            else
                legend(handles.axes3,'hide')
            end
            
        else
            %             if handles.Case == 11 %Oesophagus Case1
            %                 struct_name_case = {'GTV', 'CTVB'};
            %             end
            %
            %             if handles.Case == 31 %Prostate Case1
            %
            %             end
            
            % Retrieving information about the structures to show
            struct_name_case = handles.FIELDRTGSCases(handles.Case(1,1)).StructOverunder(handles.Case(1,2)).Cases;
            
            %% Plotting GS contour
            if get(handles.refvolume_checkbox, 'Value') == 1 % plot GS if reference checkbox is activated
                
                struct_name_temp = struct_name_case{index_selected};
                % Retrieving information about the selected structure
                GSsstructure  = getStructNum([struct_name_temp '_GS'], handles.planC, indexS);
                
                % GS's structure
                if isempty(handles.planC{handles.structfieldnum}(GSsstructure).contour)
                    xCoords = [];
                    yCoords = [];
                    %         set(handles.axes3.Children(1), 'XData', xCoords, 'YData', yCoords, 'Color', [0 1 1], 'linewidth', 1,...
                    %             'visible', 'off')
                    hhl = hhl;
                else
                    hL1 = [];
                    numSegs = length(handles.planC{handles.structfieldnum}(GSsstructure).contour(handles.currentSlice_ax).segments);
                    
                    %     thisax = handles.axes3;
                    %     plot2handle = get(thisax, 'UserData');
                    
                    yCoordsGS = [];
                    xCoordsSG = [];
                    for indexseg = 1 : numSegs
                        pointsM = handles.planC{handles.structfieldnum}(GSsstructure).contour(handles.currentSlice_ax).segments(indexseg).points;
                        if isempty(pointsM) % Where there are no contours to plot
                            xCoords = [];
                            yCoords = [];
                            %                 set(handles.axes3.Children(1), 'XData', [], 'YData', [], 'Color', [0 1 1], 'linewidth', 1,...
                            %                     'visible', 'off')
                            hhl = hhl;
                            
                        else
                            xCoords = pointsM(:,1);
                            yCoords = pointsM(:,2);
                            
                            [rowV, colV] = aapmtom(xCoords/scale, yCoords/scale, xCTOffset/scale, yCTOffset/scale, imageWidth);
                            
                            rowV = [rowV; rowV(1)];
                            colV = [colV; colV(1)];
                            
                            %Update contour plot
                            hold on
                            % hPlot = plot([colV; colV(1)], [rowV; rowV(1)], 'rs-', 'Color', [0 1 1], 'LineStyle', '-');
                            %                 set(handles.axes3.Children(1), 'XData', colV, 'YData', rowV, 'Color', [0 1 1], 'linewidth', 1,...
                            %                     'visible', 'on')
                            hL1 = plot(colV, rowV, 'Color', [1 1 0], 'linewidth', round(get(handles.Slider_zoom, 'Value')));
                            refreshdata(handles.axes3) %, 'Position', [0.1, 0.012])
                            
                            %                     imshow(blue);
                            %                 blue = cat(3, zeros(handles.volsize(1 : 2)), zeros(handles.volsize(1 : 2)), ones(handles.volsize(1 : 2)));
                            %                 set(handles.hShow.Parent.Children(7), 'AlphaData', blue);
                            %                 https://blogs.mathworks.com/steve/2009/02/18/image-overlay-using-transparency/
                            
                            %     % GS's structure
                            %     pointsGS = handles.planC{handles.structfieldnum}(4*(index_selected - 1) + 13).contour(handles.currentSlice_ax).segments(indexseg).points;
                            %     yCoordsGS = [yCoordsGS; pointsGS(:,2); pointsGS(1,2)]; % + yT; Rotation component (??)  /Users/concetta/OneDrive - Cardiff
                            %     % University/Concetta_work_cardiff/Projects/Other_projects/CERR-master/CERR_core/Viewers
                            %     % showStructures line 227
                            %     xCoordsSG = [xCoordsSG; pointsGS(:,1); pointsGS(1,1)]; %  + xT; Rotation component (??)  /Users/concetta/OneDrive - Cardiff
                            %     % University/Concetta_work_cardiff/Projects/Other_projects/CERR-master/CERR_core/Viewers
                            %     % showStructures line 227
                        end
                    end
                    
                    % Legend for the contour
                    if isempty(hL1)
                        
                    else
                        hhl = [hhl hL1];
                        labelhhl = [labelhhl {'Reference volume'}];
                    end
                    
                end
            end
            
            %% Plotting user's structure
            struct_name_temp = struct_name_case{index_selected};
            % Retrieving information about the selected structure
            USsstructure  = getStructNum(struct_name_temp, handles.planC, indexS);
            
            
            if isempty(handles.planC{handles.structfieldnum}(USsstructure).contour)
                xCoords = [];
                yCoords = [];
                %         set(handles.axes3.Children(2), 'XData', xCoords, 'YData', yCoords, 'Color', [1 0 0], 'linewidth', 1,...
                %             'visible', 'off')
                hhl = hhl;
            else
                hL2 = [];
                numSegs = length(handles.planC{handles.structfieldnum}(USsstructure).contour(handles.currentSlice_ax).segments);
                
                yCoordsGS = [];
                xCoordsSG = [];
                rowV = [];
                colV = [];
                for indexseg = 1 : numSegs
                    pointsM = handles.planC{handles.structfieldnum}(USsstructure).contour(handles.currentSlice_ax).segments(indexseg).points;
                    
                    if isempty(pointsM) % Where there are no contours to plot
                        xCoords = [];
                        yCoords = [];
                        %                 set(handles.axes3.Children(2), 'XData', xCoords, 'YData', yCoords, 'Color', [1 0 0], 'linewidth', 1,...
                        %                     'visible', 'off')
                        hhl = hhl;
                        
                    else
                        xCoords = pointsM(:,1);
                        yCoords = pointsM(:,2);
                        
                        [rowV, colV] = aapmtom(xCoords/scale, yCoords/scale, xCTOffset/scale, yCTOffset/scale, imageWidth);
                        
                        rowV = [rowV; rowV(1)];
                        colV = [colV; colV(1)];
                        
                        %Update contour plot
                        hold on
                        
                        % hPlot = plot([colV; colV(1)], [rowV; rowV(1)], 'rs-', 'Color', [1 0 0], 'LineStyle', '-');
                        %                 set(handles.axes3.Children(2), 'XData', colV, 'YData', rowV, 'Color', [1 0 0], 'linewidth', 1,...
                        %                     'visible', 'on')
                        hL2 = plot(colV, rowV, 'Color', [0 1 1], 'linewidth', round(get(handles.Slider_zoom, 'Value')));
                        refreshdata(handles.axes3)
                        
                        
                        %     % GS's structure
                        %     pointsGS = handles.planC{handles.structfieldnum}(4*(index_selected - 1) + 13).contour(handles.currentSlice_ax).segments(indexseg).points;
                        %     yCoordsGS = [yCoordsGS; pointsGS(:,2); pointsGS(1,2)]; % + yT; Rotation component (??)  /Users/concetta/OneDrive - Cardiff
                        %     % University/Concetta_work_cardiff/Projects/Other_projects/CERR-master/CERR_core/Viewers
                        %     % showStructures line 227
                        %     xCoordsSG = [xCoordsSG; pointsGS(:,1); pointsGS(1,1)]; %  + xT; Rotation component (??)  /Users/concetta/OneDrive - Cardiff
                        %     % University/Concetta_work_cardiff/Projects/Other_projects/CERR-master/CERR_core/Viewers
                        %     % showStructures line 227
                    end
                end
                
                if isempty(hL2)
                    
                else
                    hhl = [hhl, hL2];
                    labelhhl = [labelhhl, {'User volume'}];
                end
            end
            
            
            
            %% Plotting OARs
            % Case = 1; % To be changed based on the analysed case
            indexS = handles.planC{end};
            % To be changed based on the analysed case
            %             if handles.Case == 11 % Oesophagus Case1
            %                 struct_name_case_OARs = {'Vertebra_GS'; 'Aorta_GS'; 'Right lung_GS'; 'Pericardium/great vessels_GS'; ...
            %                     'Liver_GS'; 'Stomach_GS'; 'Azygous vein_GS'; 'Left main bronchus_GS'; 'Left lung_GS'};
            %             end
            %
            %             if handles.Case == 31 % Prostate Case1
            %                 struct_name_case_OARs = {'Bladder_GS'; 'Bowel_GS'; 'Lt FemHead_GS'; 'Penile bulb_GS'; 'Rectum_GS'; 'Rt FemHead_GS'};
            %             end
            
            struct_name_case_OARs = handles.FIELDRTGSCases(handles.Case(1,1)).OARs(handles.Case(1,2)).Cases;            
            colorlab = handles.FIELDRTGSCases(handles.Case(1,1)).OARs(handles.Case(1,2)).colorlab;
            
            for indexOARs = 1 : size(struct_name_case_OARs, 2) % Retrieving information about the selected structure
                OARs_struct_name_temp = [struct_name_case_OARs{indexOARs}];
                OARsstructure  = getStructNum([OARs_struct_name_temp '_GS'], handles.planC, indexS);
                
                % Structures obtained with the intersection of GS and user's structures
                if isempty(handles.planC{handles.structfieldnum}(OARsstructure).contour)
                    xCoords = [];
                    yCoords = [];
                    %         set(handles.axes3.Children(3), 'XData', xCoords, 'YData', yCoords, 'Color', [0 1 0], 'linewidth', 1,...
                    %             'visible', 'off')
                    hhl = hhl;
                else
                    hL4 = [];
                    numSegs = length(handles.planC{handles.structfieldnum}(OARsstructure).contour(handles.currentSlice_ax).segments);
                    
                    yCoordsGS = [];
                    xCoordsSG = [];
                    rowV = [];
                    colV = [];
                    for indexseg = 1 : numSegs
                        pointsM = handles.planC{handles.structfieldnum}(OARsstructure).contour(handles.currentSlice_ax).segments(indexseg).points;
                        
                        if isempty(pointsM) % Where there are no contours to plot
                            xCoords = [];
                            yCoords = [];
                            %                 set(handles.axes3.Children(3), 'XData', xCoords, 'YData', yCoords, 'Color', [0 1 0], 'linewidth', 6,...
                            %                     'visible', 'off')
                            hhl = hhl;
                            
                        else
                            xCoords = pointsM(:,1);
                            yCoords = pointsM(:,2);
                            
                            [rowV, colV] = aapmtom(xCoords/scale, yCoords/scale, xCTOffset/scale, yCTOffset/scale, imageWidth);
                            
                            rowV = [rowV; rowV(1)];
                            colV = [colV; colV(1)];
                            
                            %Update contour plot
                            hold on
                            %             hPlot = plot([colV; colV(1)], [rowV; rowV(1)], 'rs-', 'Color', [0 1 0], 'LineStyle', '-');
                            %                 set(handles.axes3.Children(3), 'XData', colV, 'YData', rowV, 'Color', [0 1 0], 'linewidth', 1,...
                            %                     'visible', 'on')
                            
                            
                            
                            hL4 = plot(colV, rowV, 'color', colorlab(indexOARs, :), 'linewidth', round(get(handles.Slider_zoom, 'Value')));
                            
                            %     % GS's structure
                            %     pointsGS = handles.planC{handles.structfieldnum}(4*(index_selected - 1) + 13).contour(handles.currentSlice_ax).segments(indexseg).points;
                            %     yCoordsGS = [yCoordsGS; pointsGS(:,2); pointsGS(1,2)]; % + yT; Rotation component (??)  /Users/concetta/OneDrive - Cardiff
                            %     % University/Concetta_work_cardiff/Projects/Other_projects/CERR-master/CERR_core/Viewers
                            %     % showStructures line 227
                            %     xCoordsSG = [xCoordsSG; pointsGS(:,1); pointsGS(1,1)]; %  + xT; Rotation component (??)  /Users/concetta/OneDrive - Cardiff
                            %     % University/Concetta_work_cardiff/Projects/Other_projects/CERR-master/CERR_core/Viewers
                            %     % showStructures line 227
                            refreshdata(handles.axes3)
                            hold on
                        end
                    end
                    
                    % Legend for the contour
                    % Only for the first segment. We don't need a label for each segment
                    if isempty(hL4)
                        
                    else
                        hhl = [hhl, hL4];
                        % struct_name_temp = struct_name_case_OARs{indexOARs};
                        % labelhhl = [labelhhl {struct_name_temp(1 : end - 3)}];
                        labelhhl = [labelhhl {struct_name_case_OARs{indexOARs}}];
                    end
                end
            end
                        
            % Checking if there is an overlapping between the user selected
            % structure and each OAR in the current slice. If yes, a flag appears.
            clear CC
            % S = [];
            %         new_handle = copyobj(hhllegend, handles.output);
            %         delete(get(new_handle,'Children'))
            if isempty(index_selected)
                
            else
                indexlegendline = 0; % flag to avoid adding two blank lines in the legend more than one time
                
                for indexoarmask = 1 : size(handles.OARsdummy_finalmask, 2) % size(handles.OARs_finalmask, 2)
                    tempmask_dummmyOARs = handles.OARsdummy_finalmask{index_selected, indexoarmask};
                    flag_label = 0; % flag to avoid adding in the legend two centroids derived from the same OAR
                    
                    % Checking if the US outline is intersection with the
                    % dummy OAR
                    CC_dummy(index_selected, indexoarmask) = bwconncomp(tempmask_dummmyOARs(:, :, handles.currentSlice_ax)); % connectivity = 8 (default)
                    if CC_dummy(index_selected, indexoarmask).NumObjects > 0
                        % there is an intersection so the centroid is
                        % computed from the intersection between the US
                        % outline and the OAR (not the dummy OAR)
                        % tempmask_OARs = handles.OARs_finalmask{index_selected, indexoarmask};
                        % CC(index_selected, indexoarmask) = bwconncomp(tempmask_OARs(:, :, handles.currentSlice_ax)); % connectivity = 8 (default)
                        % TT = regionprops(CC(index_selected, indexoarmask), 'Centroid', 'Area');
                        TT = regionprops(CC_dummy(index_selected, indexoarmask), 'Centroid', 'Area');
                        % S = [S; TT];
                        % find the maximum area
                        [val_TT, idx_TT] = max([TT.Area]);
                        
                        % To add two blank lines in the legend
                        if indexlegendline == 0
                            for indexx1 = 1 : 2
                                hL5 = plot([0.1 0.1], 'Color', 'none', 'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 1);
                                hhl = [hhl, hL5];
                                struct_name_temp =  ' ';
                                labelhhl = [labelhhl struct_name_temp];
                            end
                            indexlegendline = 1;
                        end
                        
                        hLmask_new = plot(TT(idx_TT, 1).Centroid(1, 1), TT(idx_TT, 1).Centroid(1, 2), 'LineStyle', 'none', 'Marker', '.', 'color', colorlab(indexoarmask, :), 'MarkerSize', (round(get(handles.Slider_zoom, 'Value')) * 10));
                        
                        if flag_label == 0 % so to avoid adding in the legend two centroids derived from the same OAR
                            hhl = [hhl, hLmask_new];
                            struct_name_temp = [struct_name_case_OARs{indexoarmask}];
                            labelhhl = [labelhhl , {['You have included too much ' struct_name_temp ' in your ' struct_name_case{index_selected}]}];
                            flag_label = 1;
                        end
                                
%                         for indexNumObjects = 1 : CC_dummy(index_selected, indexoarmask).NumObjects % CC(index_selected, indexoarmask).NumObjects
%                             if TT(indexNumObjects, 1).Area > 30 && flag_label == 0 % We don't want to annotate just a single point
%                                 % hAnnotAxes = findall(gcf,'type','annotation');
%                                 % arrow_coordinates_x = [(S(1, 1).Centroid(1, 1) - 60), TT(indexNumObjects, 1).Centroid(1, 1)];
%                                 % arrow_coordinates_y = [(S(1, 1).Centroid(1, 2) - 60), TT(indexNumObjects, 1).Centroid(1, 2)];
%                                 % Plots the centroid of the intersection between the US selected structure and the hidden OARs
%                                 
%                                 % To add two blank lines in the legend
%                                 if indexlegendline == 0
%                                     for indexx1 = 1 : 2
%                                         hL5 = plot([0.1 0.1], 'Color', 'none', 'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 1);
%                                         hhl = [hhl, hL5];
%                                         struct_name_temp =  ' ';
%                                         labelhhl = [labelhhl struct_name_temp];
%                                     end
%                                     indexlegendline = 1;
%                                 end
%                                 
%                                 hLmask_new = plot(TT(indexNumObjects, 1).Centroid(1, 1), TT(indexNumObjects, 1).Centroid(1, 2), 'LineStyle', 'none', 'Marker', '.', 'color', handles.colorlab(indexoarmask, :),'MarkerSize', (round(get(handles.Slider_zoom, 'Value')) * 10));
%                                 % from data space to normalized figure coordinates
%                                 % [XFIG, YFIG, DEEP] = ds2fig(arrow_coordinates_x, arrow_coordinates_y);
%                                 % if ~isempty(hAnnotAxes)
%                                 %h_annotation = annotation(handles.output, 'textarrow', XFIG, YFIG);
%                                 % else
%                                 % h_annotation = annotation(handles.output, 'textarrow', XFIG, YFIG, 'String', ['You have included .... in your GTV/CTVB/PTV']);
%                                 % end
%                                 % h_annotation.Color = 'red';
%                                 % h_annotation.FontSize = 14;
%                                 
%                                 if flag_label == 0 % so to avoid adding in the legend two centroids derived from the same OAR
%                                     hhl = [hhl, hLmask_new];
%                                     struct_name_temp = [struct_name_case_OARs{indexoarmask}];
%                                     labelhhl = [labelhhl , {['You have included too much ' struct_name_temp(1 : end - 3) ' in your ' struct_name_case{index_selected}]}];
%                                     flag_label = 1;
%                                 end
%                             else
%                                 
%                             end
%                             % text(round(S(indexNumObjects, 1).Centroid(1, 1)), round(S(indexNumObjects, 1).Centroid(1, 2)),'\rightarrow ciao polla')
%                         end
                    else
                        % Nothing appears
                    end
                end
                
                
                
                if size(hhl, 2) >= 1
                    hhllegend = legend(hhl, labelhhl, 'FontSize', 12, 'TextColor', 'white', 'Color', 'Black');
                else
                    legend(handles.axes3,'hide')
                end
                
                % a2 = axes;
                
                % %             if size(hhlmask, 2) >= 1
                % %                 a = axes('position', get(gca,'position'), 'visible', 'on');
                % %                 hhllegendmask = legend(a, hLmask_new, labelhhlmask, 'FontSize', 12, 'TextColor', 'white', 'Color', 'Black', 'Location', 'NorthWest');
                % %                 legend update
                % %             else
                % %                 a = axes('position',get(gca,'position'),'visible','on');
                % %                 % legend(a,'hide')
                % %                 legend update
                % %             end
            end
        end
        
    case 3
        if isempty(index_selected) % When no structure is selected
            legend(handles.axes3,'hide')
        else
            
            % Retrieving information about the structures to show
            struct_name_case = handles.FIELDRTGSCases(handles.Case(1,1)).StructMinmax(handles.Case(1,2)).Cases;
            GSsstructure_original  = getStructNum([struct_name_case{1, index_selected} '_GS'], handles.planC, indexS);
            USsstructure_original  = getStructNum(struct_name_case{1, index_selected}, handles.planC, indexS);
            GSsstructuremax  = getStructNum([struct_name_case{1, index_selected} '_max_GS'], handles.planC, indexS);
            GSsstructuremin  = getStructNum([struct_name_case{1, index_selected} '_min_GS'], handles.planC, indexS);
            
            %% Plotting GS structure
            if get(handles.refvolume_checkbox, 'Value') == 1 % plot GS if reference checkbox is activated
                % GS's structure
                if isempty(handles.planC{handles.structfieldnum}(USsstructure_original).contour)
                    xCoords = [];
                    yCoords = [];
                    %         set(handles.axes3.Children(1), 'XData', xCoords, 'YData', yCoords, 'Color', [0 1 1], 'linewidth', 1,...
                    %             'visible', 'off')
                    hhl = hhl;
                else
                    hL1 = [];
                    numSegs = length(handles.planC{handles.structfieldnum}(GSsstructure_original).contour(handles.currentSlice_ax).segments);
                    
                    %     thisax = handles.axes3;
                    %     plot2handle = get(thisax, 'UserData');
                    
                    yCoordsGS = [];
                    xCoordsSG = [];
                    for indexseg = 1 : numSegs
                        pointsM = handles.planC{handles.structfieldnum}(GSsstructure_original).contour(handles.currentSlice_ax).segments(indexseg).points;
                        if isempty(pointsM) % Where there are no contours to plot
                            xCoords = [];
                            yCoords = [];
                            %                 set(handles.axes3.Children(1), 'XData', [], 'YData', [], 'Color', [0 1 1], 'linewidth', 1,...
                            %                     'visible', 'off')
                            hhl = hhl;
                            
                        else
                            xCoords = pointsM(:,1);
                            yCoords = pointsM(:,2);
                            
                            [rowV, colV] = aapmtom(xCoords/scale, yCoords/scale, xCTOffset/scale, yCTOffset/scale, imageWidth);
                            
                            rowV = [rowV; rowV(1)];
                            colV = [colV; colV(1)];
                            
                            %Update contour plot
                            hold on
                            % hPlot = plot([colV; colV(1)], [rowV; rowV(1)], 'rs-', 'Color', [0 1 1], 'LineStyle', '-');
                            %                 set(handles.axes3.Children(1), 'XData', colV, 'YData', rowV, 'Color', [0 1 1], 'linewidth', 1,...
                            %                     'visible', 'on')
                            hold on
                            hL1 = plot(colV, rowV, 'Color', [1 1 0], 'linewidth', round(get(handles.Slider_zoom, 'Value')));
                            
                            refreshdata(handles.axes3) %, 'Position', [0.1, 0.012])
                            
                            %                     imshow(blue);
                            %                 blue = cat(3, zeros(handles.volsize(1 : 2)), zeros(handles.volsize(1 : 2)), ones(handles.volsize(1 : 2)));
                            %                 set(handles.hShow.Parent.Children(7), 'AlphaData', blue);
                            %                 https://blogs.mathworks.com/steve/2009/02/18/image-overlay-using-transparency/
                            
                            %     % GS's structure
                            %     pointsGS = handles.planC{handles.structfieldnum}(4*(index_selected - 1) + 13).contour(handles.currentSlice_ax).segments(indexseg).points;
                            %     yCoordsGS = [yCoordsGS; pointsGS(:,2); pointsGS(1,2)]; % + yT; Rotation component (??)  /Users/concetta/OneDrive - Cardiff
                            %     % University/Concetta_work_cardiff/Projects/Other_projects/CERR-master/CERR_core/Viewers
                            %     % showStructures line 227
                            %     xCoordsSG = [xCoordsSG; pointsGS(:,1); pointsGS(1,1)]; %  + xT; Rotation component (??)  /Users/concetta/OneDrive - Cardiff
                            %     % University/Concetta_work_cardiff/Projects/Other_projects/CERR-master/CERR_core/Viewers
                            %     % showStructures line 227
                        end
                    end
                    
                    % Legend for the contour
                    if isempty(hL1)
                        
                    else
                        hhl = [hhl hL1];
                        labelhhl = [labelhhl {'Reference volume'}];
                    end
                    
                end
            end
            
            %% Plotting user's structure
            struct_name_temp = struct_name_case{index_selected};
            % Retrieving information about the selected structure
            USsstructure  = getStructNum(struct_name_temp, handles.planC, indexS);
            
            
            if isempty(handles.planC{handles.structfieldnum}(USsstructure).contour)
                xCoords = [];
                yCoords = [];
                %         set(handles.axes3.Children(2), 'XData', xCoords, 'YData', yCoords, 'Color', [1 0 0], 'linewidth', 1,...
                %             'visible', 'off')
                hhl = hhl;
            else
                hL2 = [];
                numSegs = length(handles.planC{handles.structfieldnum}(USsstructure).contour(handles.currentSlice_ax).segments);
                
                yCoordsGS = [];
                xCoordsSG = [];
                rowV = [];
                colV = [];
                for indexseg = 1 : numSegs
                    pointsM = handles.planC{handles.structfieldnum}(USsstructure).contour(handles.currentSlice_ax).segments(indexseg).points;
                    
                    if isempty(pointsM) % Where there are no contours to plot
                        xCoords = [];
                        yCoords = [];
                        %                 set(handles.axes3.Children(2), 'XData', xCoords, 'YData', yCoords, 'Color', [1 0 0], 'linewidth', 1,...
                        %                     'visible', 'off')
                        hhl = hhl;
                        
                    else
                        xCoords = pointsM(:,1);
                        yCoords = pointsM(:,2);
                        
                        [rowV, colV] = aapmtom(xCoords/scale, yCoords/scale, xCTOffset/scale, yCTOffset/scale, imageWidth);
                        
                        rowV = [rowV; rowV(1)];
                        colV = [colV; colV(1)];
                        
                        %Update contour plot
                        hold on
                        
                        % hPlot = plot([colV; colV(1)], [rowV; rowV(1)], 'rs-', 'Color', [1 0 0], 'LineStyle', '-');
                        %                 set(handles.axes3.Children(2), 'XData', colV, 'YData', rowV, 'Color', [1 0 0], 'linewidth', 1,...
                        %                     'visible', 'on')
                        hL2 = plot(colV, rowV, 'Color', [0 1 1], 'linewidth', round(get(handles.Slider_zoom, 'Value')));
                        refreshdata(handles.axes3)
                        
                        %     % GS's structure
                        %     pointsGS = handles.planC{handles.structfieldnum}(4*(index_selected - 1) + 13).contour(handles.currentSlice_ax).segments(indexseg).points;
                        %     yCoordsGS = [yCoordsGS; pointsGS(:,2); pointsGS(1,2)]; % + yT; Rotation component (??)  /Users/concetta/OneDrive - Cardiff
                        %     % University/Concetta_work_cardiff/Projects/Other_projects/CERR-master/CERR_core/Viewers
                        %     % showStructures line 227
                        %     xCoordsSG = [xCoordsSG; pointsGS(:,1); pointsGS(1,1)]; %  + xT; Rotation component (??)  /Users/concetta/OneDrive - Cardiff
                        %     % University/Concetta_work_cardiff/Projects/Other_projects/CERR-master/CERR_core/Viewers
                        %     % showStructures line 227
                    end
                end
                
                % Legend for the contour
                if isempty(hL2)
                    
                else
                    hhl = [hhl hL2];
                    labelhhl = [labelhhl {'User volume'}];
                end
            end
            
            %             % Structures obtained with the intersection of GS and user's structures
            %             if isempty(handles.planC{handles.structfieldnum}(Intersectedstructures).contour)
            %                 xCoords = [];
            %                 yCoords = [];
            %                 %         set(handles.axes3.Children(3), 'XData', xCoords, 'YData', yCoords, 'Color', [0 1 0], 'linewidth', 1,...
            %                 %             'visible', 'off')
            %                 hhl = hhl;
            %             else
            %                 hL3 = [];
            %                 numSegs = length(handles.planC{handles.structfieldnum}(Intersectedstructures).contour(handles.currentSlice_ax).segments);
            %
            %                 yCoordsGS = [];
            %                 xCoordsSG = [];
            %                 rowV = [];
            %                 colV = [];
            %                 for indexseg = 1 : numSegs
            %                     pointsM = handles.planC{handles.structfieldnum}(Intersectedstructures).contour(handles.currentSlice_ax).segments(indexseg).points;
            %
            %                     if isempty(pointsM) % Where there are no contours to plot
            %                         xCoords = [];
            %                         yCoords = [];
            %                         %                 set(handles.axes3.Children(3), 'XData', xCoords, 'YData', yCoords, 'Color', [0 1 0], 'linewidth', 6,...
            %                         %                     'visible', 'off')
            %                         hhl = hhl;
            %
            %                     else
            %                         xCoords = pointsM(:,1);
            %                         yCoords = pointsM(:,2);
            %
            %                         [rowV, colV] = aapmtom(xCoords/scale, yCoords/scale, xCTOffset/scale, yCTOffset/scale, imageWidth);
            %
            %                         rowV = [rowV; rowV(1)];
            %                         colV = [colV; colV(1)];
            %
            %                         %Update contour plot
            %                         hold on
            %                         %             hPlot = plot([colV; colV(1)], [rowV; rowV(1)], 'rs-', 'Color', [0 1 0], 'LineStyle', '-');
            %                         %                 set(handles.axes3.Children(3), 'XData', colV, 'YData', rowV, 'Color', [0 1 0], 'linewidth', 1,...
            %                         %                     'visible', 'on')
            %                         hL3 = plot(colV, rowV, 'Color', [0 1 0], 'linewidth', 1);
            %
            %                         %     % GS's structure
            %                         %     pointsGS = handles.planC{handles.structfieldnum}(4*(index_selected - 1) + 13).contour(handles.currentSlice_ax).segments(indexseg).points;
            %                         %     yCoordsGS = [yCoordsGS; pointsGS(:,2); pointsGS(1,2)]; % + yT; Rotation component (??)  /Users/concetta/OneDrive - Cardiff
            %                         %     % University/Concetta_work_cardiff/Projects/Other_projects/CERR-master/CERR_core/Viewers
            %                         %     % showStructures line 227
            %                         %     xCoordsSG = [xCoordsSG; pointsGS(:,1); pointsGS(1,1)]; %  + xT; Rotation component (??)  /Users/concetta/OneDrive - Cardiff
            %                         %     % University/Concetta_work_cardiff/Projects/Other_projects/CERR-master/CERR_core/Viewers
            %                         %     % showStructures line 227
            %                         refreshdata(handles.axes3)
            %                         hold on
            %                     end
            %                 end
            %
            %                 % Legend for the contour
            %                 if isempty(hL3)
            %
            %                 else
            %                     hhl = [hhl hL3];
            %                     labelhhl = [labelhhl {'Common region'}];
            %                 end
            %             end
            
            % if get(handles.refvolume_checkbox, 'Value') == 1 % plot GS's structure max if reference checkbox is activated
            
            %% Plotting GS's max structure
            if isempty(handles.planC{handles.structfieldnum}(GSsstructuremax).contour)
                xCoords = [];
                yCoords = [];
                %         set(handles.axes3.Children(3), 'XData', xCoords, 'YData', yCoords, 'Color', [0 1 0], 'linewidth', 1,...
                %             'visible', 'off')
                hhl = hhl;
            else
                hL4 = [];
                numSegs = length(handles.planC{handles.structfieldnum}(GSsstructuremax).contour(handles.currentSlice_ax).segments);
                
                yCoordsGS = [];
                xCoordsSG = [];
                rowV = [];
                colV = [];
                for indexseg = 1 : numSegs
                    pointsM = handles.planC{handles.structfieldnum}(GSsstructuremax).contour(handles.currentSlice_ax).segments(indexseg).points;
                    
                    if isempty(pointsM) % Where there are no contours to plot
                        xCoords = [];
                        yCoords = [];
                        %                 set(handles.axes3.Children(3), 'XData', xCoords, 'YData', yCoords, 'Color', [0 1 0], 'linewidth', 6,...
                        %                     'visible', 'off')
                        hhl = hhl;
                        
                    else
                        xCoords = pointsM(:,1);
                        yCoords = pointsM(:,2);
                        
                        [rowV, colV] = aapmtom(xCoords/scale, yCoords/scale, xCTOffset/scale, yCTOffset/scale, imageWidth);
                        
                        rowV = [rowV; rowV(1)];
                        colV = [colV; colV(1)];
                        
                        %Update contour plot
                        hold on
                        %             hPlot = plot([colV; colV(1)], [rowV; rowV(1)], 'rs-', 'Color', [0 1 0], 'LineStyle', '-');
                        %                 set(handles.axes3.Children(3), 'XData', colV, 'YData', rowV, 'Color', [0 1 0], 'linewidth', 1,...
                        %                     'visible', 'on')
                        % hL4 = plot(colV, rowV, 'Color', [0 1 1], 'linewidth', 1, 'LineStyle', '--');
                        hL4 = plot(colV, rowV, 'Color', [1 0 0], 'linewidth', round(get(handles.Slider_zoom, 'Value')));
                        
                        %     % GS's structure
                        %     pointsGS = handles.planC{handles.structfieldnum}(4*(index_selected - 1) + 13).contour(handles.currentSlice_ax).segments(indexseg).points;
                        %     yCoordsGS = [yCoordsGS; pointsGS(:,2); pointsGS(1,2)]; % + yT; Rotation component (??)  /Users/concetta/OneDrive - Cardiff
                        %     % University/Concetta_work_cardiff/Projects/Other_projects/CERR-master/CERR_core/Viewers
                        %     % showStructures line 227
                        %     xCoordsSG = [xCoordsSG; pointsGS(:,1); pointsGS(1,1)]; %  + xT; Rotation component (??)  /Users/concetta/OneDrive - Cardiff
                        %     % University/Concetta_work_cardiff/Projects/Other_projects/CERR-master/CERR_core/Viewers
                        %     % showStructures line 227
                        refreshdata(handles.axes3)
                        hold on
                    end
                end
                
                % Legend for the contour
                if isempty(hL4)
                    
                else
                    hhl = [hhl hL4];
                    struct_name_temp = struct_name_case{index_selected};
                    labelhhl = [labelhhl {[struct_name_temp ' max']}];
                end
            end
            %  end
            
            
            % if get(handles.refvolume_checkbox, 'Value') == 1 % plot GS's structure min if reference checkbox is activated
            
            %% Plotting GS's min structure
            if isempty(handles.planC{handles.structfieldnum}(GSsstructuremin).contour)
                xCoords = [];
                yCoords = [];
                %         set(handles.axes3.Children(3), 'XData', xCoords, 'YData', yCoords, 'Color', [0 1 0], 'linewidth', 1,...
                %             'visible', 'off')
                hhl = hhl;
            else
                hL5 = [];
                numSegs = length(handles.planC{handles.structfieldnum}(GSsstructuremin).contour(handles.currentSlice_ax).segments);
                
                yCoordsGS = [];
                xCoordsSG = [];
                rowV = [];
                colV = [];
                for indexseg = 1 : numSegs
                    pointsM = handles.planC{handles.structfieldnum}(GSsstructuremin).contour(handles.currentSlice_ax).segments(indexseg).points;
                    
                    if isempty(pointsM) % Where there are no contours to plot
                        xCoords = [];
                        yCoords = [];
                        %                 set(handles.axes3.Children(3), 'XData', xCoords, 'YData', yCoords, 'Color', [0 1 0], 'linewidth', 6,...
                        %                     'visible', 'off')
                        hhl = hhl;
                        
                    else
                        xCoords = pointsM(:,1);
                        yCoords = pointsM(:,2);
                        
                        [rowV, colV] = aapmtom(xCoords/scale, yCoords/scale, xCTOffset/scale, yCTOffset/scale, imageWidth);
                        
                        rowV = [rowV; rowV(1)];
                        colV = [colV; colV(1)];
                        
                        %Update contour plot
                        hold on
                        %             hPlot = plot([colV; colV(1)], [rowV; rowV(1)], 'rs-', 'Color', [0 1 0], 'LineStyle', '-');
                        %                 set(handles.axes3.Children(3), 'XData', colV, 'YData', rowV, 'Color', [0 1 0], 'linewidth', 1,...
                        %                     'visible', 'on')
                        % hL5 = plot(colV, rowV, 'Color', [1 1 0], 'linewidth', 2, 'LineStyle', ':');
                        hL5 = plot(colV, rowV, 'Color', [0 1 0], 'linewidth', round(get(handles.Slider_zoom, 'Value')));
                        
                        %     % GS's structure
                        %     pointsGS = handles.planC{handles.structfieldnum}(4*(index_selected - 1) + 13).contour(handles.currentSlice_ax).segments(indexseg).points;
                        %     yCoordsGS = [yCoordsGS; pointsGS(:,2); pointsGS(1,2)]; % + yT; Rotation component (??)  /Users/concetta/OneDrive - Cardiff
                        %     % University/Concetta_work_cardiff/Projects/Other_projects/CERR-master/CERR_core/Viewers
                        %     % showStructures line 227
                        %     xCoordsSG = [xCoordsSG; pointsGS(:,1); pointsGS(1,1)]; %  + xT; Rotation component (??)  /Users/concetta/OneDrive - Cardiff
                        %     % University/Concetta_work_cardiff/Projects/Other_projects/CERR-master/CERR_core/Viewers
                        %     % showStructures line 227
                        refreshdata(handles.axes3)
                        hold on
                    end
                end
                
                % Legend for the contour
                if isempty(hL5)
                    
                else
                    hhl = [hhl hL5];
                    struct_name_temp = struct_name_case{index_selected};
                    labelhhl = [labelhhl {[struct_name_temp ' min']}];
                end
            end
            % end
            
            %             if size(hhl, 2) == 3
            %                 hh = legend(hhl, {'Under region', 'Over region', 'Common region'}, 'FontSize', 12, 'TextColor', 'white', 'Color', 'Black');
            %             end
            %             if size(hhl, 2) == 2
            %                 hh = legend(hhl, {'Under region', 'Over region'}, 'FontSize', 12, 'TextColor', 'white', 'Color', 'Black');
            %             end
            %             if size(hhl, 2) == 1
            %                 if labelhhl(1, 1) == 1
            %                     hh = legend(hhl, {'Under region'}, 'FontSize', 12, 'TextColor', 'white', 'Color', 'Black');
            %                 else
            %                     hh = legend(hhl, {'Over region'}, 'FontSize', 12, 'TextColor', 'white', 'Color', 'Black');
            %                 end
            %             end
            %             if size(hhl, 2) < 1
            %                 legend(handles.axes3,'hide')
            %             end
            if size(hhl, 2) >= 1
                hhllegend = legend(hhl, labelhhl, 'FontSize', 12, 'TextColor', 'white', 'Color', 'Black');
            else
                legend(handles.axes3,'hide')
            end
        end
        
    case 4
        if isempty(index_selected) % When no structure is selected
            legend(handles.axes3,'hide')
        else
            
            if get(handles.refvolume_checkbox, 'Value') == 1 % plot GS if reference checkbox is activated
                % Retrieving information about the structures to show
                struct_name_case = handles.FIELDRTGSCases(handles.Case(1,1)).Structures(handles.Case(1,2)).Cases;
                GSsstructure_original  = getStructNum([struct_name_case{1, index_selected} '_GS'], handles.planC, indexS);
                
                %% Plotting GS structure
                if isempty(handles.planC{handles.structfieldnum}(GSsstructure_original).contour)
                    xCoords = [];
                    yCoords = [];
                    %         set(handles.axes3.Children(1), 'XData', xCoords, 'YData', yCoords, 'Color', [0 1 1], 'linewidth', 1,...
                    %             'visible', 'off')
                    hhl = hhl;
                else
                    hL1 = [];
                    numSegs = length(handles.planC{handles.structfieldnum}(GSsstructure_original).contour(handles.currentSlice_ax).segments);
                    
                    %     thisax = handles.axes3;
                    %     plot2handle = get(thisax, 'UserData');
                    
                    yCoordsGS = [];
                    xCoordsSG = [];
                    for indexseg = 1 : numSegs
                        pointsM = handles.planC{handles.structfieldnum}(GSsstructure_original).contour(handles.currentSlice_ax).segments(indexseg).points;
                        if isempty(pointsM) % Where there are no contours to plot
                            xCoords = [];
                            yCoords = [];
                            %                 set(handles.axes3.Children(1), 'XData', [], 'YData', [], 'Color', [0 1 1], 'linewidth', 1,...
                            %                     'visible', 'off')
                            hhl = hhl;
                            
                        else
                            xCoords = pointsM(:,1);
                            yCoords = pointsM(:,2);
                            
                            [rowV, colV] = aapmtom(xCoords/scale, yCoords/scale, xCTOffset/scale, yCTOffset/scale, imageWidth);
                            
                            rowV = [rowV; rowV(1)];
                            colV = [colV; colV(1)];
                            
                            %Update contour plot
                            hold on
                            % hPlot = plot([colV; colV(1)], [rowV; rowV(1)], 'rs-', 'Color', [0 1 1], 'LineStyle', '-');
                            %                 set(handles.axes3.Children(1), 'XData', colV, 'YData', rowV, 'Color', [0 1 1], 'linewidth', 1,...
                            %                     'visible', 'on')
                            hL1 = plot(colV, rowV, 'Color', [1 1 0], 'linewidth', round(get(handles.Slider_zoom, 'Value')));
                            
                            refreshdata(handles.axes3)                            
                            hold on
                            
                            %                     imshow(blue);
                            %                 blue = cat(3, zeros(handles.volsize(1 : 2)), zeros(handles.volsize(1 : 2)), ones(handles.volsize(1 : 2)));
                            %                 set(handles.hShow.Parent.Children(7), 'AlphaData', blue);
                            %                 https://blogs.mathworks.com/steve/2009/02/18/image-overlay-using-transparency/
                            
                            %     % GS's structure
                            %     pointsGS = handles.planC{handles.structfieldnum}(4*(index_selected - 1) + 13).contour(handles.currentSlice_ax).segments(indexseg).points;
                            %     yCoordsGS = [yCoordsGS; pointsGS(:,2); pointsGS(1,2)]; % + yT; Rotation component (??)  /Users/concetta/OneDrive - Cardiff
                            %     % University/Concetta_work_cardiff/Projects/Other_projects/CERR-master/CERR_core/Viewers
                            %     % showStructures line 227
                            %     xCoordsSG = [xCoordsSG; pointsGS(:,1); pointsGS(1,1)]; %  + xT; Rotation component (??)  /Users/concetta/OneDrive - Cardiff
                            %     % University/Concetta_work_cardiff/Projects/Other_projects/CERR-master/CERR_core/Viewers
                            %     % showStructures line 227
                        end
                    end
                    
                    % Legend for the contour
                    if isempty(hL1)
                        
                    else
                        hhl = [hhl hL1];
                        labelhhl = [labelhhl {'Reference volume'}];
                    end
                end
            end
            
            %% Plotting user's structure
            
            % Retrieving information about the structures to show
            struct_name_case = handles.FIELDRTGSCases(handles.Case(1,1)).Structures(handles.Case(1,2)).Cases;
            USsstructure_original  = getStructNum(struct_name_case{1, index_selected}, handles.planC, indexS);
            
            if isempty(handles.planC{handles.structfieldnum}(USsstructure_original).contour)
                xCoords = [];
                yCoords = [];
                %         set(handles.axes3.Children(2), 'XData', xCoords, 'YData', yCoords, 'Color', [1 0 0], 'linewidth', 1,...
                %             'visible', 'off')
                hhl = hhl;
            else
                hL2 = [];
                numSegs = length(handles.planC{handles.structfieldnum}(USsstructure_original).contour(handles.currentSlice_ax).segments);
                
                yCoordsGS = [];
                xCoordsSG = [];
                rowV = [];
                colV = [];
                for indexseg = 1 : numSegs
                    pointsM = handles.planC{handles.structfieldnum}(USsstructure_original).contour(handles.currentSlice_ax).segments(indexseg).points;
                    
                    if isempty(pointsM) % Where there are no contours to plot
                        xCoords = [];
                        yCoords = [];
                        %                 set(handles.axes3.Children(2), 'XData', xCoords, 'YData', yCoords, 'Color', [1 0 0], 'linewidth', 1,...
                        %                     'visible', 'off')
                        hhl = hhl;
                        
                    else
                        xCoords = pointsM(:,1);
                        yCoords = pointsM(:,2);
                        
                        [rowV, colV] = aapmtom(xCoords/scale, yCoords/scale, xCTOffset/scale, yCTOffset/scale, imageWidth);
                        
                        rowV = [rowV; rowV(1)];
                        colV = [colV; colV(1)];
                        
                        %Update contour plot
                        hold on
                        
                        % hPlot = plot([colV; colV(1)], [rowV; rowV(1)], 'rs-', 'Color', [1 0 0], 'LineStyle', '-');
                        %                 set(handles.axes3.Children(2), 'XData', colV, 'YData', rowV, 'Color', [1 0 0], 'linewidth', 1,...
                        %                     'visible', 'on')
                        hL2 = plot(colV, rowV, 'Color', [0 1 1], 'linewidth', round(get(handles.Slider_zoom, 'Value')));
                        refreshdata(handles.axes3)
                        
                        %     % GS's structure
                        %     pointsGS = handles.planC{handles.structfieldnum}(4*(index_selected - 1) + 13).contour(handles.currentSlice_ax).segments(indexseg).points;
                        %     yCoordsGS = [yCoordsGS; pointsGS(:,2); pointsGS(1,2)]; % + yT; Rotation component (??)  /Users/concetta/OneDrive - Cardiff
                        %     % University/Concetta_work_cardiff/Projects/Other_projects/CERR-master/CERR_core/Viewers
                        %     % showStructures line 227
                        %     xCoordsSG = [xCoordsSG; pointsGS(:,1); pointsGS(1,1)]; %  + xT; Rotation component (??)  /Users/concetta/OneDrive - Cardiff
                        %     % University/Concetta_work_cardiff/Projects/Other_projects/CERR-master/CERR_core/Viewers
                        %     % showStructures line 227
                    end
                end
                
                % Legend for the contour
                if isempty(hL2)
                    
                else
                    hhl = [hhl hL2];
                    labelhhl = [labelhhl {'User volume'}];
                end
            end
            
            if size(hhl, 2) >= 1
                hhllegend = legend(hhl, labelhhl, 'FontSize', 12, 'TextColor', 'white', 'Color', 'Black');
            else
                legend(handles.axes3,'hide')
            end
        end
    case 5
        if isempty(index_selected) % When no structure is selected
            legend(handles.axes3,'hide')
        else
            %             % Legend for the contour
            %             if isempty(hL1)
            %
            %             else
            %             end
        end
end

% %% Sagittal view
%
% % Setting the axes where visualize the contours
% axes(handles.axes5)
%
% switch n
%
%     case 1
%         if isempty(index_selected) % When no structure is selected
%             legend(handles.axes5,'hide')
%         else
%             % GS's structure
%
%             % The contour will be extracted from the sagittal/coronal mask and plotted
%             %
%             sag_mask_bound = [];
%             sag_mask_bound = bwboundaries(handles.sagmasks_structs(GSsstructure).mask(:, :, handles.currentSlice_sag));
%
%
%             Xrange = handles.sagittalxlim(2) - handles.sagittalxlim(1);
%             Yrange = handles.sagittalylim(2) - handles.sagittalylim(1);
%
%
%             if isempty(sag_mask_bound)
%                 xCoords = [];
%                 yCoords = [];
%                 %         set(handles.axes5.Children(1), 'XData', xCoords, 'YData', yCoords, 'Color', [0 1 1], 'linewidth', 1,...
%                 %             'visible', 'off')
%                 hhl = hhl;
%             else
%                 hL1 = [];
%                 numSegs = length(sag_mask_bound);
%
%                 %     thisax = handles.axes5;
%                 %     plot2handle = get(thisax, 'UserData');
%
%                 yCoordsGS = [];
%                 xCoordsSG = [];
%                 for indexseg = 1 : numSegs
%
%                     X = []; Y = [];
%                     X = (sag_mask_bound{1, 1}(:, 1) * (handles.sagittalxlim(2) - handles.sagittalxlim(1))) / size(handles.sagmasks_structs(GSsstructure).mask, 1) + handles.sagittalxlim(1);
%                     Y = (sag_mask_bound{1, 1}(:, 2) * (handles.sagittalxlim(2) - handles.sagittalxlim(1))) / size(handles.sagmasks_structs(GSsstructure).mask, 2) + handles.sagittalxlim(1);
%
%                     pointsM = [X Y];
%
%                     if isempty(pointsM) % Where there are no contours to plot
%                         xCoords = [];
%                         yCoords = [];
%                         %                 set(handles.axes5.Children(1), 'XData', [], 'YData', [], 'Color', [0 1 1], 'linewidth', 1,...
%                         %                     'visible', 'off')
%                         hhl = hhl;
%
%                     else
%
%                         rowV = [pointsM(:, 1); pointsM(1, 1)];
%                         colV = [pointsM(:, 2); pointsM(1, 2)];
%
%                         %Update contour plot
%                         hold on
%                         % hPlot = plot([colV; colV(1)], [rowV; rowV(1)], 'rs-', 'Color', [0 1 1], 'LineStyle', '-');
%                         %                 set(handles.axes5.Children(1), 'XData', colV, 'YData', rowV, 'Color', [0 1 1], 'linewidth', 1,...
%                         %                     'visible', 'on')
%                         hL1 = plot(colV, rowV, 'Color', [1 1 0], 'linewidth', round(get(handles.Slider_zoom, 'Value')));
%                         refreshdata(handles.axes5) %, 'Position', [0.1, 0.012])
%
%                         %                     imshow(blue);
%                         %                 blue = cat(3, zeros(handles.volsize(1 : 2)), zeros(handles.volsize(1 : 2)), ones(handles.volsize(1 : 2)));
%                         %                 set(handles.hShow.Parent.Children(7), 'AlphaData', blue);
%                         %                 https://blogs.mathworks.com/steve/2009/02/18/image-overlay-using-transparency/
%
%                         %     % GS's structure
%                         %     pointsGS = handles.planC{handles.structfieldnum}(4*(index_selected - 1) + 13).contour(handles.currentSlice_ax).segments(indexseg).points;
%                         %     yCoordsGS = [yCoordsGS; pointsGS(:,2); pointsGS(1,2)]; % + yT; Rotation component (??)  /Users/concetta/OneDrive - Cardiff
%                         %     % University/Concetta_work_cardiff/Projects/Other_projects/CERR-master/CERR_core/Viewers
%                         %     % showStructures line 227
%                         %     xCoordsSG = [xCoordsSG; pointsGS(:,1); pointsGS(1,1)]; %  + xT; Rotation component (??)  /Users/concetta/OneDrive - Cardiff
%                         %     % University/Concetta_work_cardiff/Projects/Other_projects/CERR-master/CERR_core/Viewers
%                         %     % showStructures line 227
%                     end
%                 end
%
% %                 % Legend for the contour
% %                 if isempty(hL1)
% %
% %                 else
% %                     hhl = [hhl hL1];
% %                     labelhhl = [labelhhl {'Under region'}];
% %                 end
%
%             end
%
%             % User's structure
%
%             % The contour will be extracted from the sagittal/coronal mask and plotted
%             %
%             sag_mask_bound = [];
%             sag_mask_bound = bwboundaries(handles.sagmasks_structs(USsstructure).mask(:, :, handles.currentSlice_sag));
%
%
%             Xrange = handles.sagittalxlim(2) - handles.sagittalxlim(1);
%             Yrange = handles.sagittalylim(2) - handles.sagittalylim(1);
%
%             if isempty(sag_mask_bound)
%                 xCoords = [];
%                 yCoords = [];
%                 %         set(handles.axes5.Children(2), 'XData', xCoords, 'YData', yCoords, 'Color', [1 0 0], 'linewidth', 1,...
%                 %             'visible', 'off')
%                 hhl = hhl;
%             else
%                 hL2 = [];
%                 numSegs = length(sag_mask_bound);
%
%                 yCoordsGS = [];
%                 xCoordsSG = [];
%                 rowV = [];
%                 colV = [];
%                 for indexseg = 1 : numSegs
%
%                     X = []; Y = [];
%                     X = (sag_mask_bound{1, 1}(:, 1) * (handles.sagittalxlim(2) - handles.sagittalxlim(1))) / size(handles.sagmasks_structs(USsstructure).mask, 1) + handles.sagittalxlim(1);
%                     Y = (sag_mask_bound{1, 1}(:, 2) * (handles.sagittalxlim(2) - handles.sagittalxlim(1))) / size(handles.sagmasks_structs(USsstructure).mask, 2) + handles.sagittalxlim(1);
%
%                     pointsM = [X Y];
%
%                     if isempty(pointsM) % Where there are no contours to plot
%                         xCoords = [];
%                         yCoords = [];
%                         %                 set(handles.axes5.Children(2), 'XData', xCoords, 'YData', yCoords, 'Color', [1 0 0], 'linewidth', 1,...
%                         %                     'visible', 'off')
%                         hhl = hhl;
%
%                     else
%                         rowV = [pointsM(:, 1); pointsM(1, 1)];
%                         colV = [pointsM(:, 2); pointsM(1, 2)];
%
%                         %Update contour plot
%                         hold on
%
%                         % hPlot = plot([colV; colV(1)], [rowV; rowV(1)], 'rs-', 'Color', [1 0 0], 'LineStyle', '-');
%                         %                 set(handles.axes5.Children(2), 'XData', colV, 'YData', rowV, 'Color', [1 0 0], 'linewidth', 1,...
%                         %                     'visible', 'on')
%
%                         hL2 = plot(colV, rowV, 'Color', [0 1 1], 'linewidth', round(get(handles.Slider_zoom, 'Value')));
%                         refreshdata(handles.axes5)
%
%                         %     % GS's structure
%                         %     pointsGS = handles.planC{handles.structfieldnum}(4*(index_selected - 1) + 13).contour(handles.currentSlice_ax).segments(indexseg).points;
%                         %     yCoordsGS = [yCoordsGS; pointsGS(:,2); pointsGS(1,2)]; % + yT; Rotation component (??)  /Users/concetta/OneDrive - Cardiff
%                         %     % University/Concetta_work_cardiff/Projects/Other_projects/CERR-master/CERR_core/Viewers
%                         %     % showStructures line 227
%                         %     xCoordsSG = [xCoordsSG; pointsGS(:,1); pointsGS(1,1)]; %  + xT; Rotation component (??)  /Users/concetta/OneDrive - Cardiff
%                         %     % University/Concetta_work_cardiff/Projects/Other_projects/CERR-master/CERR_core/Viewers
%                         %     % showStructures line 227
%                     end
%                 end
%
% %                 % Legend for the contour
% %                 if isempty(hL2)
% %
% %                 else
% %                     hhl = [hhl hL2];
% %                     labelhhl = [labelhhl {'Over region'}];
% %                 end
%             end
%
%             % Structures obtained with the intersection of GS and user's structures
%
%             % The contour will be extracted from the sagittal/coronal mask and plotted
%             %
%             sag_mask_bound = [];
%             sag_mask_bound = bwboundaries(handles.sagmasks_structs(Intersectedstructures).mask(:, :, handles.currentSlice_sag));
%
%
%             Xrange = handles.sagittalxlim(2) - handles.sagittalxlim(1);
%             Yrange = handles.sagittalylim(2) - handles.sagittalylim(1);
%
%             if isempty(sag_mask_bound)
%                 xCoords = [];
%                 yCoords = [];
%                 %         set(handles.axes5.Children(3), 'XData', xCoords, 'YData', yCoords, 'Color', [0 1 0], 'linewidth', 1,...
%                 %             'visible', 'off')
%                 hhl = hhl;
%             else
%                 hL3 = [];
%                 numSegs = length(sag_mask_bound);
%
%                 yCoordsGS = [];
%                 xCoordsSG = [];
%                 rowV = [];
%                 colV = [];
%                 for indexseg = 1 : numSegs
%
%                     X = []; Y = [];
%                     X = (sag_mask_bound{1, 1}(:, 1) * (handles.sagittalxlim(2) - handles.sagittalxlim(1))) / size(handles.sagmasks_structs(Intersectedstructures).mask, 1) + handles.sagittalxlim(1);
%                     Y = (sag_mask_bound{1, 1}(:, 2) * (handles.sagittalxlim(2) - handles.sagittalxlim(1))) / size(handles.sagmasks_structs(Intersectedstructures).mask, 2) + handles.sagittalxlim(1);
%
%                     pointsM = [X Y];
%
%                     if isempty(pointsM) % Where there are no contours to plot
%                         xCoords = [];
%                         yCoords = [];
%                         %                 set(handles.axes5.Children(3), 'XData', xCoords, 'YData', yCoords, 'Color', [0 1 0], 'linewidth', 6,...
%                         %                     'visible', 'off')
%                         hhl = hhl;
%
%                     else
%                         rowV = [pointsM(:, 1); pointsM(1, 1)];
%                         colV = [pointsM(:, 2); pointsM(1, 2)];
%
%                         %Update contour plot
%                         hold on
%                         %             hPlot = plot([colV; colV(1)], [rowV; rowV(1)], 'rs-', 'Color', [0 1 0], 'LineStyle', '-');
%                         %                 set(handles.axes5.Children(3), 'XData', colV, 'YData', rowV, 'Color', [0 1 0], 'linewidth', 1,...
%                         %                     'visible', 'on')
%                         hL3 = plot(colV, rowV, 'Color', [0 1 0], 'linewidth', round(get(handles.Slider_zoom, 'Value')));
%
%                         %     % GS's structure
%                         %     pointsGS = handles.planC{handles.structfieldnum}(4*(index_selected - 1) + 13).contour(handles.currentSlice_ax).segments(indexseg).points;
%                         %     yCoordsGS = [yCoordsGS; pointsGS(:,2); pointsGS(1,2)]; % + yT; Rotation component (??)  /Users/concetta/OneDrive - Cardiff
%                         %     % University/Concetta_work_cardiff/Projects/Other_projects/CERR-master/CERR_core/Viewers
%                         %     % showStructures line 227
%                         %     xCoordsSG = [xCoordsSG; pointsGS(:,1); pointsGS(1,1)]; %  + xT; Rotation component (??)  /Users/concetta/OneDrive - Cardiff
%                         %     % University/Concetta_work_cardiff/Projects/Other_projects/CERR-master/CERR_core/Viewers
%                         %     % showStructures line 227
%                         refreshdata(handles.axes5)
%                         hold on
%                     end
%                 end
%
% %                 % Legend for the contour
% %                 if isempty(hL3)
% %
% %                 else
% %                     hhl = [hhl hL3];
% %                     labelhhl = [labelhhl {'Common region'}];
% %                 end
%             end
%
%             %             if size(hhl, 2) == 3
%             %                 hh = legend(hhl, {'Under region', 'Over region', 'Common region'}, 'FontSize', 12, 'TextColor', 'white', 'Color', 'Black');
%             %             end
%             %             if size(hhl, 2) == 2
%             %                 hh = legend(hhl, {'Under region', 'Over region'}, 'FontSize', 12, 'TextColor', 'white', 'Color', 'Black');
%             %             end
%             %             if size(hhl, 2) == 1
%             %                 if labelhhl(1, 1) == 1
%             %                     hh = legend(hhl, {'Under region'}, 'FontSize', 12, 'TextColor', 'white', 'Color', 'Black');
%             %                 else
%             %                     hh = legend(hhl, {'Over region'}, 'FontSize', 12, 'TextColor', 'white', 'Color', 'Black');
%             %                 end
%             %             end
%             %             if size(hhl, 2) < 1
%             %                 legend(handles.axes5,'hide')
%             %             end
% %             if size(hhl, 2) >= 1
% %                 hhllegend = legend(hhl, labelhhl, 'FontSize', 12, 'TextColor', 'white', 'Color', 'Black');
% %             else
% %                 legend(handles.axes5,'hide')
% %             end
%         end
%
%     case 2
%         if isempty(index_selected) % When no structure is selected
%             legend(handles.axes5,'hide')
%         else
%             struct_name_case = {'GTV', 'CTVB'}; %, 'PTV'};
%             % User's structure
%             struct_name_temp = struct_name_case{index_selected};
%             % Retrieving information about the selected structure
%             USsstructure  = getStructNum(struct_name_temp, handles.planC, indexS);
%
%             % The contour will be extracted from the sagittal/coronal mask and plotted
%             %
%             sag_mask_bound = [];
%             sag_mask_bound = bwboundaries(handles.sagmasks_structs(USsstructure).mask(:, :, handles.currentSlice_sag));
%
%
%             Xrange = handles.sagittalxlim(2) - handles.sagittalxlim(1);
%             Yrange = handles.sagittalylim(2) - handles.sagittalylim(1);
%
%             if isempty(sag_mask_bound)
%                 xCoords = [];
%                 yCoords = [];
%                 %         set(handles.axes5.Children(2), 'XData', xCoords, 'YData', yCoords, 'Color', [1 0 0], 'linewidth', 1,...
%                 %             'visible', 'off')
%                 hhl = hhl;
%             else
%                 hL2 = [];
%                 numSegs = length(sag_mask_bound);
%
%                 yCoordsGS = [];
%                 xCoordsSG = [];
%                 rowV = [];
%                 colV = [];
%                 for indexseg = 1 : numSegs
%
%                     X = []; Y = [];
%                     X = (sag_mask_bound{1, 1}(:, 1) * (handles.sagittalxlim(2) - handles.sagittalxlim(1))) / size(handles.sagmasks_structs(USsstructure).mask, 1) + handles.sagittalxlim(1);
%                     Y = (sag_mask_bound{1, 1}(:, 2) * (handles.sagittalxlim(2) - handles.sagittalxlim(1))) / size(handles.sagmasks_structs(USsstructure).mask, 2) + handles.sagittalxlim(1);
%
%                     pointsM = [X Y];
%
%                     if isempty(pointsM) % Where there are no contours to plot
%                         xCoords = [];
%                         yCoords = [];
%                         %                 set(handles.axes5.Children(2), 'XData', xCoords, 'YData', yCoords, 'Color', [1 0 0], 'linewidth', 1,...
%                         %                     'visible', 'off')
%                         hhl = hhl;
%
%                     else
%                         rowV = [pointsM(:, 1); pointsM(1, 1)];
%                         colV = [pointsM(:, 2); pointsM(1, 2)];
%
%                         %Update contour plot
%                         hold on
%
%                         % hPlot = plot([colV; colV(1)], [rowV; rowV(1)], 'rs-', 'Color', [1 0 0], 'LineStyle', '-');
%                         %                 set(handles.axes5.Children(2), 'XData', colV, 'YData', rowV, 'Color', [1 0 0], 'linewidth', 1,...
%                         %                     'visible', 'on')
%                         hL2 = plot(colV, rowV, 'Color', [0 1 1], 'linewidth', round(get(handles.Slider_zoom, 'Value')));
%                         refreshdata(handles.axes5)
%
%
%                         %     % GS's structure
%                         %     pointsGS = handles.planC{handles.structfieldnum}(4*(index_selected - 1) + 13).contour(handles.currentSlice_ax).segments(indexseg).points;
%                         %     yCoordsGS = [yCoordsGS; pointsGS(:,2); pointsGS(1,2)]; % + yT; Rotation component (??)  /Users/concetta/OneDrive - Cardiff
%                         %     % University/Concetta_work_cardiff/Projects/Other_projects/CERR-master/CERR_core/Viewers
%                         %     % showStructures line 227
%                         %     xCoordsSG = [xCoordsSG; pointsGS(:,1); pointsGS(1,1)]; %  + xT; Rotation component (??)  /Users/concetta/OneDrive - Cardiff
%                         %     % University/Concetta_work_cardiff/Projects/Other_projects/CERR-master/CERR_core/Viewers
%                         %     % showStructures line 227
%                     end
%                 end
%
% %                 if isempty(hL2)
% %
% %                 else
% %                     hhl = [hhl, hL2];
% %                     labelhhl = [labelhhl, {'User volume'}];
% %                 end
%             end
%
%
%         end
%
%         Case = 1; % To be changed based on the analysed case
%         indexS = handles.planC{end};
%         % To be changed based on the analysed case
%         if Case == 1
%             struct_name_case_OARs = {'Vertebra_GS'; 'Aorta_GS'; 'Right lung_GS'; 'Pericardium/great vessels_GS'; ...
%                 'Liver_GS'; 'Stomach_GS'; 'Azygous vein_GS'; 'Left main bronchus_GS'; 'Left lung_GS'};
%
%         end
%         for indexOARs = 1 : size(struct_name_case_OARs, 1) % Retrieving information about the selected structure
%             struct_name_temp = [struct_name_case_OARs{indexOARs}];
%             OARsstructure  = getStructNum(struct_name_temp, handles.planC, indexS);
%
%             % The contour will be extracted from the sagittal/coronal mask and plotted
%             %
%             sag_mask_bound = [];
%             sag_mask_bound = bwboundaries(handles.sagmasks_structs(OARsstructure).mask(:, :, handles.currentSlice_sag));
%
%
%             Xrange = handles.sagittalxlim(2) - handles.sagittalxlim(1);
%             Yrange = handles.sagittalylim(2) - handles.sagittalylim(1);
%
%             if isempty(sag_mask_bound)
%                 xCoords = [];
%                 yCoords = [];
%                 %         set(handles.axes5.Children(3), 'XData', xCoords, 'YData', yCoords, 'Color', [0 1 0], 'linewidth', 1,...
%                 %             'visible', 'off')
%                 hhl = hhl;
%             else
%                 hL4 = [];
%                 numSegs = length(sag_mask_bound);
%
%                 yCoordsGS = [];
%                 xCoordsSG = [];
%                 rowV = [];
%                 colV = [];
%                 for indexseg = 1 : numSegs
%
%                     X = []; Y = [];
%                     X = (sag_mask_bound{1, 1}(:, 1) * (handles.sagittalxlim(2) - handles.sagittalxlim(1))) / size(handles.sagmasks_structs(OARsstructure).mask, 1) + handles.sagittalxlim(1);
%                     Y = (sag_mask_bound{1, 1}(:, 2) * (handles.sagittalxlim(2) - handles.sagittalxlim(1))) / size(handles.sagmasks_structs(OARsstructure).mask, 2) + handles.sagittalxlim(1);
%
%                     pointsM = [X Y];
%
%                     if isempty(pointsM) % Where there are no contours to plot
%                         xCoords = [];
%                         yCoords = [];
%                         %                 set(handles.axes5.Children(3), 'XData', xCoords, 'YData', yCoords, 'Color', [0 1 0], 'linewidth', 6,...
%                         %                     'visible', 'off')
%                         hhl = hhl;
%
%                     else
%
%                         rowV = [pointsM(:, 1); pointsM(1, 1)];
%                         colV = [pointsM(:, 2); pointsM(1, 2)];
%
%                         rowV = [rowV; rowV(1)];
%                         colV = [colV; colV(1)];
%
%                         %Update contour plot
%                         hold on
%                         %             hPlot = plot([colV; colV(1)], [rowV; rowV(1)], 'rs-', 'Color', [0 1 0], 'LineStyle', '-');
%                         %                 set(handles.axes5.Children(3), 'XData', colV, 'YData', rowV, 'Color', [0 1 0], 'linewidth', 1,...
%                         %                     'visible', 'on')
%
%
%
%                         hL4 = plot(colV, rowV, 'color', handles.colorlab(indexOARs, :), 'linewidth', round(get(handles.Slider_zoom, 'Value')));
%
%                         %     % GS's structure
%                         %     pointsGS = handles.planC{handles.structfieldnum}(4*(index_selected - 1) + 13).contour(handles.currentSlice_ax).segments(indexseg).points;
%                         %     yCoordsGS = [yCoordsGS; pointsGS(:,2); pointsGS(1,2)]; % + yT; Rotation component (??)  /Users/concetta/OneDrive - Cardiff
%                         %     % University/Concetta_work_cardiff/Projects/Other_projects/CERR-master/CERR_core/Viewers
%                         %     % showStructures line 227
%                         %     xCoordsSG = [xCoordsSG; pointsGS(:,1); pointsGS(1,1)]; %  + xT; Rotation component (??)  /Users/concetta/OneDrive - Cardiff
%                         %     % University/Concetta_work_cardiff/Projects/Other_projects/CERR-master/CERR_core/Viewers
%                         %     % showStructures line 227
%                         refreshdata(handles.axes5)
%                         hold on
%                     end
%                 end
%
%                 % Legend for the contour
%                 % Only for the first segment. We don't need a label for each segment
% %                 if isempty(hL4)
% %
% %                 else
% %                     hhl = [hhl, hL4];
% %                     struct_name_temp = struct_name_case_OARs{indexOARs};
% %                     labelhhl = [labelhhl {struct_name_temp(1 : end - 3)}];
% %                 end
%             end
%         end
%
% %         % To add two blank lines in the legend
% %         for indexx1 = 1 : 2
% %             hL5 = plot([0.1 0.1], [0.1 0.1],'k','LineWidth', 0.1);
% %             hhl = [hhl, hL5];
% %             struct_name_temp =  ' ';
% %             labelhhl = [labelhhl struct_name_temp];
% %         end
% %
% %         % Checking if there is an overlapping between the user selected
% %         % structure and each OAR in the current slice. If yes, a flag appears.
% %         clear CC
% %         S = [];
% %         %         new_handle = copyobj(hhllegend, handles.output);
%         %         delete(get(new_handle,'Children'))
%         if isempty(index_selected)
%
%         else
% %             for indexoarmask = 1 : size(handles.OARs_finalmask, 2)
% %                 tempmask = handles.OARs_finalmask{index_selected, indexoarmask};
% %                 flag_label = 0; % flag to avoid adding in the legend two centroids derived from the same OAR
% %                 % Checking if there is an overlapping
% %                 CC(index_selected, indexoarmask) = bwconncomp(tempmask(:, :, handles.currentSlice_ax)); % connectivity = 8 (default)
% %                 if CC(index_selected, indexoarmask).NumObjects > 0
% %                     TT = regionprops(CC(index_selected, indexoarmask), 'Centroid', 'Area');
% %                     S = [S; TT];
% %                     for indexNumObjects = 1 : CC(index_selected, indexoarmask).NumObjects
% %                         if TT(indexNumObjects, 1).Area > 30 % We don't want to annotate just a single point
% %                             % hAnnotAxes = findall(gcf,'type','annotation');
% %                             % arrow_coordinates_x = [(S(1, 1).Centroid(1, 1) - 60), TT(indexNumObjects, 1).Centroid(1, 1)];
% %                             % arrow_coordinates_y = [(S(1, 1).Centroid(1, 2) - 60), TT(indexNumObjects, 1).Centroid(1, 2)];
% %                             % Plots the centroid of the intersection between the US selected structure and the hidden OARs
% %                             hLmask_new = plot(TT(indexNumObjects, 1).Centroid(1, 1), TT(indexNumObjects, 1).Centroid(1, 2), 'Marker', '.', 'color', handles.colorlab(indexoarmask, :),'MarkerSize', 20);
% %                             % from data space to normalized figure coordinates
% %                             % [XFIG, YFIG, DEEP] = ds2fig(arrow_coordinates_x, arrow_coordinates_y);
% %                             % if ~isempty(hAnnotAxes)
% %                             %h_annotation = annotation(handles.output, 'textarrow', XFIG, YFIG);
% %                             % else
% %                             % h_annotation = annotation(handles.output, 'textarrow', XFIG, YFIG, 'String', ['You have included .... in your GTV/CTVB/PTV']);
% %                             % end
% %                             % h_annotation.Color = 'red';
% %                             % h_annotation.FontSize = 14;
% %
% %                             if flag_label == 0 % so to avoid adding in the legend two centroids derived from the same OAR
% %                                 hhl = [hhl, hLmask_new];
% %                                 struct_name_temp = [struct_name_case_OARs{indexoarmask}];
% %                                 labelhhl = [labelhhl , {['You have included too much ' struct_name_temp(1 : end - 3) ' in your ' struct_name_case{index_selected}]}];
% %                                 flag_label = 1;
% %                             end
% %                         else
% %
% %                         end
% %                         % text(round(S(indexNumObjects, 1).Centroid(1, 1)), round(S(indexNumObjects, 1).Centroid(1, 2)),'\rightarrow ciao polla')
% %                     end
% %                 else
% %                     % Nothing appears
% %                 end
% %             end
%
%
% %
% %             if size(hhl, 2) >= 1
% %                 hhllegend = legend(hhl, labelhhl, 'FontSize', 12, 'TextColor', 'white', 'Color', 'Black');
% %             else
% %                 legend(handles.axes5,'hide')
% %             end
%
%             % a2 = axes;
%
%             % %             if size(hhlmask, 2) >= 1
%             % %                 a = axes('position', get(gca,'position'), 'visible', 'on');
%             % %                 hhllegendmask = legend(a, hLmask_new, labelhhlmask, 'FontSize', 12, 'TextColor', 'white', 'Color', 'Black', 'Location', 'NorthWest');
%             % %                 legend update
%             % %             else
%             % %                 a = axes('position',get(gca,'position'),'visible','on');
%             % %                 % legend(a,'hide')
%             % %                 legend update
%             % %             end
%         end
%
%     case 3
%         if isempty(index_selected) % When no structure is selected
%             legend(handles.axes5,'hide')
%         else
%             % GS's structure
%
%             % The contour will be extracted from the sagittal/coronal mask and plotted
%             %
%             sag_mask_bound = [];
%             sag_mask_bound = bwboundaries(handles.sagmasks_structs(GSsstructure_original).mask(:, :, handles.currentSlice_sag));
%
%
%             Xrange = handles.sagittalxlim(2) - handles.sagittalxlim(1);
%             Yrange = handles.sagittalylim(2) - handles.sagittalylim(1);
%
%             if isempty(sag_mask_bound)
%                 xCoords = [];
%                 yCoords = [];
%                 %         set(handles.axes5.Children(1), 'XData', xCoords, 'YData', yCoords, 'Color', [0 1 1], 'linewidth', 1,...
%                 %             'visible', 'off')
%                 hhl = hhl;
%             else
%                 hL1 = [];
%                 numSegs = length(sag_mask_bound);
%
%                 %     thisax = handles.axes5;
%                 %     plot2handle = get(thisax, 'UserData');
%
%                 yCoordsGS = [];
%                 xCoordsSG = [];
%                 for indexseg = 1 : numSegs
%
%                     X = []; Y = [];
%                     X = (sag_mask_bound{1, 1}(:, 1) * (handles.sagittalxlim(2) - handles.sagittalxlim(1))) / size(handles.sagmasks_structs(GSsstructure_original).mask, 1) + handles.sagittalxlim(1);
%                     Y = (sag_mask_bound{1, 1}(:, 2) * (handles.sagittalxlim(2) - handles.sagittalxlim(1))) / size(handles.sagmasks_structs(GSsstructure_original).mask, 2) + handles.sagittalxlim(1);
%
%                     pointsM = [X Y];
%
%                     if isempty(pointsM) % Where there are no contours to plot
%                         xCoords = [];
%                         yCoords = [];
%                         %                 set(handles.axes5.Children(1), 'XData', [], 'YData', [], 'Color', [0 1 1], 'linewidth', 1,...
%                         %                     'visible', 'off')
%                         hhl = hhl;
%
%                     else
%
%                         rowV = [pointsM(:, 1); pointsM(1, 1)];
%                         colV = [pointsM(:, 2); pointsM(1, 2)];
%
%                         rowV = [rowV; rowV(1)];
%                         colV = [colV; colV(1)];
%
%                         %Update contour plot
%                         hold on
%                         % hPlot = plot([colV; colV(1)], [rowV; rowV(1)], 'rs-', 'Color', [0 1 1], 'LineStyle', '-');
%                         %                 set(handles.axes5.Children(1), 'XData', colV, 'YData', rowV, 'Color', [0 1 1], 'linewidth', 1,...
%                         %                     'visible', 'on')
%                         hold on
%                         hL1 = plot(colV, rowV, 'Color', [1 1 0], 'linewidth', round(get(handles.Slider_zoom, 'Value')));
%
%                         refreshdata(handles.axes5) %, 'Position', [0.1, 0.012])
%
%                         %                     imshow(blue);
%                         %                 blue = cat(3, zeros(handles.volsize(1 : 2)), zeros(handles.volsize(1 : 2)), ones(handles.volsize(1 : 2)));
%                         %                 set(handles.hShow.Parent.Children(7), 'AlphaData', blue);
%                         %                 https://blogs.mathworks.com/steve/2009/02/18/image-overlay-using-transparency/
%
%                         %     % GS's structure
%                         %     pointsGS = handles.planC{handles.structfieldnum}(4*(index_selected - 1) + 13).contour(handles.currentSlice_ax).segments(indexseg).points;
%                         %     yCoordsGS = [yCoordsGS; pointsGS(:,2); pointsGS(1,2)]; % + yT; Rotation component (??)  /Users/concetta/OneDrive - Cardiff
%                         %     % University/Concetta_work_cardiff/Projects/Other_projects/CERR-master/CERR_core/Viewers
%                         %     % showStructures line 227
%                         %     xCoordsSG = [xCoordsSG; pointsGS(:,1); pointsGS(1,1)]; %  + xT; Rotation component (??)  /Users/concetta/OneDrive - Cardiff
%                         %     % University/Concetta_work_cardiff/Projects/Other_projects/CERR-master/CERR_core/Viewers
%                         %     % showStructures line 227
%                     end
%                 end
%
%                 % Legend for the contour
%                 if isempty(hL1)
%
%                 else
%                     hhl = [hhl hL1];
%                     labelhhl = [labelhhl {'Reference volume'}];
%                 end
%
%             end
%
%             % User's structure
%             struct_name_temp = struct_name_case{index_selected};
%             % Retrieving information about the selected structure
%             USsstructure  = getStructNum(struct_name_temp, handles.planC, indexS);
%
%             % The contour will be extracted from the sagittal/coronal mask and plotted
%             %
%             sag_mask_bound = [];
%             sag_mask_bound = bwboundaries(handles.sagmasks_structs(USsstructure).mask(:, :, handles.currentSlice_sag));
%
%
%             Xrange = handles.sagittalxlim(2) - handles.sagittalxlim(1);
%             Yrange = handles.sagittalylim(2) - handles.sagittalylim(1);
%
%             if isempty(sag_mask_bound)
%                 xCoords = [];
%                 yCoords = [];
%                 %         set(handles.axes5.Children(2), 'XData', xCoords, 'YData', yCoords, 'Color', [1 0 0], 'linewidth', 1,...
%                 %             'visible', 'off')
%                 hhl = hhl;
%             else
%                 hL2 = [];
%                 numSegs = length(handles.planC{handles.structfieldnum}(USsstructure).contour(handles.currentSlice_ax).segments);
%
%                 yCoordsGS = [];
%                 xCoordsSG = [];
%                 rowV = [];
%                 colV = [];
%                 for indexseg = 1 : numSegs
%
%                     X = []; Y = [];
%                     X = (sag_mask_bound{1, 1}(:, 1) * (handles.sagittalxlim(2) - handles.sagittalxlim(1))) / size(handles.sagmasks_structs(USsstructure).mask, 1) + handles.sagittalxlim(1);
%                     Y = (sag_mask_bound{1, 1}(:, 2) * (handles.sagittalxlim(2) - handles.sagittalxlim(1))) / size(handles.sagmasks_structs(USsstructure).mask, 2) + handles.sagittalxlim(1);
%
%                     pointsM = [X Y];
%
%                     if isempty(pointsM) % Where there are no contours to plot
%                         xCoords = [];
%                         yCoords = [];
%                         %                 set(handles.axes5.Children(2), 'XData', xCoords, 'YData', yCoords, 'Color', [1 0 0], 'linewidth', 1,...
%                         %                     'visible', 'off')
%                         hhl = hhl;
%
%                     else
%
%                         rowV = [pointsM(:, 1); pointsM(1, 1)];
%                         colV = [pointsM(:, 2); pointsM(1, 2)];
%
%                         rowV = [rowV; rowV(1)];
%                         colV = [colV; colV(1)];
%
%                         %Update contour plot
%                         hold on
%
%                         % hPlot = plot([colV; colV(1)], [rowV; rowV(1)], 'rs-', 'Color', [1 0 0], 'LineStyle', '-');
%                         %                 set(handles.axes5.Children(2), 'XData', colV, 'YData', rowV, 'Color', [1 0 0], 'linewidth', 1,...
%                         %                     'visible', 'on')
%                         hL2 = plot(colV, rowV, 'Color', [0 1 1], 'linewidth', round(get(handles.Slider_zoom, 'Value')));
%                         refreshdata(handles.axes5)
%
%                         %     % GS's structure
%                         %     pointsGS = handles.planC{handles.structfieldnum}(4*(index_selected - 1) + 13).contour(handles.currentSlice_ax).segments(indexseg).points;
%                         %     yCoordsGS = [yCoordsGS; pointsGS(:,2); pointsGS(1,2)]; % + yT; Rotation component (??)  /Users/concetta/OneDrive - Cardiff
%                         %     % University/Concetta_work_cardiff/Projects/Other_projects/CERR-master/CERR_core/Viewers
%                         %     % showStructures line 227
%                         %     xCoordsSG = [xCoordsSG; pointsGS(:,1); pointsGS(1,1)]; %  + xT; Rotation component (??)  /Users/concetta/OneDrive - Cardiff
%                         %     % University/Concetta_work_cardiff/Projects/Other_projects/CERR-master/CERR_core/Viewers
%                         %     % showStructures line 227
%                     end
%                 end
%
% %                 % Legend for the contour
% %                 if isempty(hL2)
% %
% %                 else
% %                     hhl = [hhl hL2];
% %                     labelhhl = [labelhhl {'User volume'}];
% %                 end
%             end
%
%             %             % Structures obtained with the intersection of GS and user's structures
%             %             if isempty(handles.planC{handles.structfieldnum}(Intersectedstructures).contour)
%             %                 xCoords = [];
%             %                 yCoords = [];
%             %                 %         set(handles.axes5.Children(3), 'XData', xCoords, 'YData', yCoords, 'Color', [0 1 0], 'linewidth', 1,...
%             %                 %             'visible', 'off')
%             %                 hhl = hhl;
%             %             else
%             %                 hL3 = [];
%             %                 numSegs = length(handles.planC{handles.structfieldnum}(Intersectedstructures).contour(handles.currentSlice_ax).segments);
%             %
%             %                 yCoordsGS = [];
%             %                 xCoordsSG = [];
%             %                 rowV = [];
%             %                 colV = [];
%             %                 for indexseg = 1 : numSegs
%             %                     pointsM = handles.planC{handles.structfieldnum}(Intersectedstructures).contour(handles.currentSlice_ax).segments(indexseg).points;
%             %
%             %                     if isempty(pointsM) % Where there are no contours to plot
%             %                         xCoords = [];
%             %                         yCoords = [];
%             %                         %                 set(handles.axes5.Children(3), 'XData', xCoords, 'YData', yCoords, 'Color', [0 1 0], 'linewidth', 6,...
%             %                         %                     'visible', 'off')
%             %                         hhl = hhl;
%             %
%             %                     else
%             %                         xCoords = pointsM(:,1);
%             %                         yCoords = pointsM(:,2);
%             %
%             %                         [rowV, colV] = aapmtom(xCoords/scale, yCoords/scale, xCTOffset/scale, yCTOffset/scale, imageWidth);
%             %
%             %                         rowV = [rowV; rowV(1)];
%             %                         colV = [colV; colV(1)];
%             %
%             %                         %Update contour plot
%             %                         hold on
%             %                         %             hPlot = plot([colV; colV(1)], [rowV; rowV(1)], 'rs-', 'Color', [0 1 0], 'LineStyle', '-');
%             %                         %                 set(handles.axes5.Children(3), 'XData', colV, 'YData', rowV, 'Color', [0 1 0], 'linewidth', 1,...
%             %                         %                     'visible', 'on')
%             %                         hL3 = plot(colV, rowV, 'Color', [0 1 0], 'linewidth', 1);
%             %
%             %                         %     % GS's structure
%             %                         %     pointsGS = handles.planC{handles.structfieldnum}(4*(index_selected - 1) + 13).contour(handles.currentSlice_ax).segments(indexseg).points;
%             %                         %     yCoordsGS = [yCoordsGS; pointsGS(:,2); pointsGS(1,2)]; % + yT; Rotation component (??)  /Users/concetta/OneDrive - Cardiff
%             %                         %     % University/Concetta_work_cardiff/Projects/Other_projects/CERR-master/CERR_core/Viewers
%             %                         %     % showStructures line 227
%             %                         %     xCoordsSG = [xCoordsSG; pointsGS(:,1); pointsGS(1,1)]; %  + xT; Rotation component (??)  /Users/concetta/OneDrive - Cardiff
%             %                         %     % University/Concetta_work_cardiff/Projects/Other_projects/CERR-master/CERR_core/Viewers
%             %                         %     % showStructures line 227
%             %                         refreshdata(handles.axes5)
%             %                         hold on
%             %                     end
%             %                 end
%             %
%             %                 % Legend for the contour
%             %                 if isempty(hL3)
%             %
%             %                 else
%             %                     hhl = [hhl hL3];
%             %                     labelhhl = [labelhhl {'Common region'}];
%             %                 end
%             %             end
%
%             % GS's structure max
%
%             % The contour will be extracted from the sagittal/coronal mask and plotted
%             %
%             sag_mask_bound = [];
%             sag_mask_bound = bwboundaries(handles.sagmasks_structs(GSsstructuremax).mask(:, :, handles.currentSlice_sag));
%
%
%             Xrange = handles.sagittalxlim(2) - handles.sagittalxlim(1);
%             Yrange = handles.sagittalylim(2) - handles.sagittalylim(1);
%
%             if isempty(sag_mask_bound)
%                 xCoords = [];
%                 yCoords = [];
%                 %         set(handles.axes5.Children(3), 'XData', xCoords, 'YData', yCoords, 'Color', [0 1 0], 'linewidth', 1,...
%                 %             'visible', 'off')
%                 hhl = hhl;
%             else
%                 hL4 = [];
%                 numSegs = length(sag_mask_bound);
%
%                 yCoordsGS = [];
%                 xCoordsSG = [];
%                 rowV = [];
%                 colV = [];
%                 for indexseg = 1 : numSegs
%
%                     X = []; Y = [];
%                     X = (sag_mask_bound{1, 1}(:, 1) * (handles.sagittalxlim(2) - handles.sagittalxlim(1))) / size(handles.sagmasks_structs(GSsstructuremax).mask, 1) + handles.sagittalxlim(1);
%                     Y = (sag_mask_bound{1, 1}(:, 2) * (handles.sagittalxlim(2) - handles.sagittalxlim(1))) / size(handles.sagmasks_structs(GSsstructuremax).mask, 2) + handles.sagittalxlim(1);
%
%                     pointsM = [X Y];
%
%                     if isempty(pointsM) % Where there are no contours to plot
%                         xCoords = [];
%                         yCoords = [];
%                         %                 set(handles.axes5.Children(3), 'XData', xCoords, 'YData', yCoords, 'Color', [0 1 0], 'linewidth', 6,...
%                         %                     'visible', 'off')
%                         hhl = hhl;
%
%                     else
%
%                         rowV = [pointsM(:, 1); pointsM(1, 1)];
%                         colV = [pointsM(:, 2); pointsM(1, 2)];
%
%                         rowV = [rowV; rowV(1)];
%                         colV = [colV; colV(1)];
%
%                         %Update contour plot
%                         hold on
%                         %             hPlot = plot([colV; colV(1)], [rowV; rowV(1)], 'rs-', 'Color', [0 1 0], 'LineStyle', '-');
%                         %                 set(handles.axes5.Children(3), 'XData', colV, 'YData', rowV, 'Color', [0 1 0], 'linewidth', 1,...
%                         %                     'visible', 'on')
%                         % hL4 = plot(colV, rowV, 'Color', [0 1 1], 'linewidth', 1, 'LineStyle', '--');
%                         hL4 = plot(colV, rowV, 'Color', [1 0 0], 'linewidth', round(get(handles.Slider_zoom, 'Value')));
%
%                         %     % GS's structure
%                         %     pointsGS = handles.planC{handles.structfieldnum}(4*(index_selected - 1) + 13).contour(handles.currentSlice_ax).segments(indexseg).points;
%                         %     yCoordsGS = [yCoordsGS; pointsGS(:,2); pointsGS(1,2)]; % + yT; Rotation component (??)  /Users/concetta/OneDrive - Cardiff
%                         %     % University/Concetta_work_cardiff/Projects/Other_projects/CERR-master/CERR_core/Viewers
%                         %     % showStructures line 227
%                         %     xCoordsSG = [xCoordsSG; pointsGS(:,1); pointsGS(1,1)]; %  + xT; Rotation component (??)  /Users/concetta/OneDrive - Cardiff
%                         %     % University/Concetta_work_cardiff/Projects/Other_projects/CERR-master/CERR_core/Viewers
%                         %     % showStructures line 227
%                         refreshdata(handles.axes5)
%                         hold on
%                     end
%                 end
%
% %                 % Legend for the contour
% %                 if isempty(hL4)
% %
% %                 else
% %                     hhl = [hhl hL4];
% %                     struct_name_temp = struct_name_case{index_selected};
% %                     labelhhl = [labelhhl {[struct_name_temp ' max']}];
% %                 end
%             end
%
%             % GS's structure min
%
%             % The contour will be extracted from the sagittal/coronal mask and plotted
%             %
%             sag_mask_bound = [];
%             sag_mask_bound = bwboundaries(handles.sagmasks_structs(GSsstructuremin).mask(:, :, handles.currentSlice_sag));
%
%
%             Xrange = handles.sagittalxlim(2) - handles.sagittalxlim(1);
%             Yrange = handles.sagittalylim(2) - handles.sagittalylim(1);
%
%             if isempty(sag_mask_bound)
%                 xCoords = [];
%                 yCoords = [];
%                 %         set(handles.axes5.Children(3), 'XData', xCoords, 'YData', yCoords, 'Color', [0 1 0], 'linewidth', 1,...
%                 %             'visible', 'off')
%                 hhl = hhl;
%             else
%                 hL5 = [];
%                 numSegs = length(handles.planC{handles.structfieldnum}(GSsstructuremin).contour(handles.currentSlice_ax).segments);
%
%                 yCoordsGS = [];
%                 xCoordsSG = [];
%                 rowV = [];
%                 colV = [];
%                 for indexseg = 1 : numSegs
%
%                     X = []; Y = [];
%                     X = (sag_mask_bound{1, 1}(:, 1) * (handles.sagittalxlim(2) - handles.sagittalxlim(1))) / size(handles.sagmasks_structs(GSsstructuremin).mask, 1) + handles.sagittalxlim(1);
%                     Y = (sag_mask_bound{1, 1}(:, 2) * (handles.sagittalxlim(2) - handles.sagittalxlim(1))) / size(handles.sagmasks_structs(GSsstructuremin).mask, 2) + handles.sagittalxlim(1);
%
%                     pointsM = [X Y];
%
%                     if isempty(pointsM) % Where there are no contours to plot
%                         xCoords = [];
%                         yCoords = [];
%                         %                 set(handles.axes5.Children(3), 'XData', xCoords, 'YData', yCoords, 'Color', [0 1 0], 'linewidth', 6,...
%                         %                     'visible', 'off')
%                         hhl = hhl;
%
%                     else
%
%                         rowV = [pointsM(:, 1); pointsM(1, 1)];
%                         colV = [pointsM(:, 2); pointsM(1, 2)];
%
%                         rowV = [rowV; rowV(1)];
%                         colV = [colV; colV(1)];
%
%                         %Update contour plot
%                         hold on
%                         %             hPlot = plot([colV; colV(1)], [rowV; rowV(1)], 'rs-', 'Color', [0 1 0], 'LineStyle', '-');
%                         %                 set(handles.axes5.Children(3), 'XData', colV, 'YData', rowV, 'Color', [0 1 0], 'linewidth', 1,...
%                         %                     'visible', 'on')
%                         % hL5 = plot(colV, rowV, 'Color', [1 1 0], 'linewidth', 2, 'LineStyle', ':');
%                         hL5 = plot(colV, rowV, 'Color', [0 1 0], 'linewidth', round(get(handles.Slider_zoom, 'Value')));
%
%                         %     % GS's structure
%                         %     pointsGS = handles.planC{handles.structfieldnum}(4*(index_selected - 1) + 13).contour(handles.currentSlice_ax).segments(indexseg).points;
%                         %     yCoordsGS = [yCoordsGS; pointsGS(:,2); pointsGS(1,2)]; % + yT; Rotation component (??)  /Users/concetta/OneDrive - Cardiff
%                         %     % University/Concetta_work_cardiff/Projects/Other_projects/CERR-master/CERR_core/Viewers
%                         %     % showStructures line 227
%                         %     xCoordsSG = [xCoordsSG; pointsGS(:,1); pointsGS(1,1)]; %  + xT; Rotation component (??)  /Users/concetta/OneDrive - Cardiff
%                         %     % University/Concetta_work_cardiff/Projects/Other_projects/CERR-master/CERR_core/Viewers
%                         %     % showStructures line 227
%                         refreshdata(handles.axes5)
%                         hold on
%                     end
%                 end
%
% %                 % Legend for the contour
% %                 if isempty(hL5)
% %
% %                 else
% %                     hhl = [hhl hL5];
% %                     struct_name_temp = struct_name_case{index_selected};
% %                     labelhhl = [labelhhl {[struct_name_temp ' min']}];
% %                 end
%             end
%
%             %             if size(hhl, 2) == 3
%             %                 hh = legend(hhl, {'Under region', 'Over region', 'Common region'}, 'FontSize', 12, 'TextColor', 'white', 'Color', 'Black');
%             %             end
%             %             if size(hhl, 2) == 2
%             %                 hh = legend(hhl, {'Under region', 'Over region'}, 'FontSize', 12, 'TextColor', 'white', 'Color', 'Black');
%             %             end
%             %             if size(hhl, 2) == 1
%             %                 if labelhhl(1, 1) == 1
%             %                     hh = legend(hhl, {'Under region'}, 'FontSize', 12, 'TextColor', 'white', 'Color', 'Black');
%             %                 else
%             %                     hh = legend(hhl, {'Over region'}, 'FontSize', 12, 'TextColor', 'white', 'Color', 'Black');
%             %                 end
%             %             end
%             %             if size(hhl, 2) < 1
%             %                 legend(handles.axes5,'hide')
%             %             end
% %             if size(hhl, 2) >= 1
% %                 hhllegend = legend(hhl, labelhhl, 'FontSize', 12, 'TextColor', 'white', 'Color', 'Black');
% %             else
% %                 legend(handles.axes5,'hide')
% %             end
%         end
%
%     case 4
%         if isempty(index_selected) % When no structure is selected
%             legend(handles.axes5,'hide')
%         else
%             % GS's structure
%
%             % The contour will be extracted from the sagittal/coronal mask and plotted
%             %
%             sag_mask_bound = [];
%             sag_mask_bound = bwboundaries(handles.sagmasks_structs(GSsstructure_original).mask(:, :, handles.currentSlice_sag));
%
%
%             Xrange = handles.sagittalxlim(2) - handles.sagittalxlim(1);
%             Yrange = handles.sagittalylim(2) - handles.sagittalylim(1);
%
%             if isempty(sag_mask_bound)
%                 xCoords = [];
%                 yCoords = [];
%                 %         set(handles.axes5.Children(1), 'XData', xCoords, 'YData', yCoords, 'Color', [0 1 1], 'linewidth', 1,...
%                 %             'visible', 'off')
%                 hhl = hhl;
%             else
%                 hL1 = [];
%                 numSegs = length(handles.planC{handles.structfieldnum}(GSsstructure_original).contour(handles.currentSlice_ax).segments);
%
%                 %     thisax = handles.axes5;
%                 %     plot2handle = get(thisax, 'UserData');
%
%                 yCoordsGS = [];
%                 xCoordsSG = [];
%                 for indexseg = 1 : numSegs
%
%                     X = []; Y = [];
%                     X = (sag_mask_bound{1, 1}(:, 1) * (handles.sagittalxlim(2) - handles.sagittalxlim(1))) / size(handles.sagmasks_structs(GSsstructure_original).mask, 1) + handles.sagittalxlim(1);
%                     Y = (sag_mask_bound{1, 1}(:, 2) * (handles.sagittalxlim(2) - handles.sagittalxlim(1))) / size(handles.sagmasks_structs(GSsstructure_original).mask, 2) + handles.sagittalxlim(1);
%
%                     pointsM = [X Y];
%
%                     if isempty(pointsM) % Where there are no contours to plot
%                         xCoords = [];
%                         yCoords = [];
%                         %                 set(handles.axes5.Children(1), 'XData', [], 'YData', [], 'Color', [0 1 1], 'linewidth', 1,...
%                         %                     'visible', 'off')
%                         hhl = hhl;
%
%                     else
%
%                         rowV = [pointsM(:, 1); pointsM(1, 1)];
%                         colV = [pointsM(:, 2); pointsM(1, 2)];
%
%                         rowV = [rowV; rowV(1)];
%                         colV = [colV; colV(1)];
%
%                         %Update contour plot
%                         hold on
%                         % hPlot = plot([colV; colV(1)], [rowV; rowV(1)], 'rs-', 'Color', [0 1 1], 'LineStyle', '-');
%                         %                 set(handles.axes5.Children(1), 'XData', colV, 'YData', rowV, 'Color', [0 1 1], 'linewidth', 1,...
%                         %                     'visible', 'on')
%                         hL1 = plot(colV, rowV, 'Color', [0 1 1], 'linewidth', round(get(handles.Slider_zoom, 'Value')));
%                         refreshdata(handles.axes5) %, 'Position', [0.1, 0.012])
%
%                         %                     imshow(blue);
%                         %                 blue = cat(3, zeros(handles.volsize(1 : 2)), zeros(handles.volsize(1 : 2)), ones(handles.volsize(1 : 2)));
%                         %                 set(handles.hShow.Parent.Children(7), 'AlphaData', blue);
%                         %                 https://blogs.mathworks.com/steve/2009/02/18/image-overlay-using-transparency/
%
%                         %     % GS's structure
%                         %     pointsGS = handles.planC{handles.structfieldnum}(4*(index_selected - 1) + 13).contour(handles.currentSlice_ax).segments(indexseg).points;
%                         %     yCoordsGS = [yCoordsGS; pointsGS(:,2); pointsGS(1,2)]; % + yT; Rotation component (??)  /Users/concetta/OneDrive - Cardiff
%                         %     % University/Concetta_work_cardiff/Projects/Other_projects/CERR-master/CERR_core/Viewers
%                         %     % showStructures line 227
%                         %     xCoordsSG = [xCoordsSG; pointsGS(:,1); pointsGS(1,1)]; %  + xT; Rotation component (??)  /Users/concetta/OneDrive - Cardiff
%                         %     % University/Concetta_work_cardiff/Projects/Other_projects/CERR-master/CERR_core/Viewers
%                         %     % showStructures line 227
%                     end
%                 end
%
% %                 % Legend for the contour
% %                 if isempty(hL1)
% %
% %                 else
% %                     hhl = [hhl hL1];
% %                     labelhhl = [labelhhl {'Reference volume'}];
% %                 end
%
%             end
%
%
%             % User's structure
%
%             % The contour will be extracted from the sagittal/coronal mask and plotted
%             %
%             sag_mask_bound = [];
%             sag_mask_bound = bwboundaries(handles.sagmasks_structs(USsstructure_original).mask(:, :, handles.currentSlice_sag));
%
%
%             Xrange = handles.sagittalxlim(2) - handles.sagittalxlim(1);
%             Yrange = handles.sagittalylim(2) - handles.sagittalylim(1);
%
%             if isempty(sag_mask_bound)
%                 xCoords = [];
%                 yCoords = [];
%                 %         set(handles.axes5.Children(2), 'XData', xCoords, 'YData', yCoords, 'Color', [1 0 0], 'linewidth', 1,...
%                 %             'visible', 'off')
%                 hhl = hhl;
%             else
%                 hL2 = [];
%                 numSegs = length(handles.planC{handles.structfieldnum}(USsstructure_original).contour(handles.currentSlice_ax).segments);
%
%                 yCoordsGS = [];
%                 xCoordsSG = [];
%                 rowV = [];
%                 colV = [];
%                 for indexseg = 1 : numSegs
%
%                     X = []; Y = [];
%                     X = (sag_mask_bound{1, 1}(:, 1) * (handles.sagittalxlim(2) - handles.sagittalxlim(1))) / size(handles.sagmasks_structs(USsstructure_original).mask, 1) + handles.sagittalxlim(1);
%                     Y = (sag_mask_bound{1, 1}(:, 2) * (handles.sagittalxlim(2) - handles.sagittalxlim(1))) / size(handles.sagmasks_structs(USsstructure_original).mask, 2) + handles.sagittalxlim(1);
%
%                     pointsM = [X Y];
%
%                     if isempty(pointsM) % Where there are no contours to plot
%                         xCoords = [];
%                         yCoords = [];
%                         %                 set(handles.axes5.Children(2), 'XData', xCoords, 'YData', yCoords, 'Color', [1 0 0], 'linewidth', 1,...
%                         %                     'visible', 'off')
%                         hhl = hhl;
%
%                     else
%
%                         rowV = [pointsM(:, 1); pointsM(1, 1)];
%                         colV = [pointsM(:, 2); pointsM(1, 2)];
%
%                         rowV = [rowV; rowV(1)];
%                         colV = [colV; colV(1)];
%
%
%                         %Update contour plot
%                         hold on
%
%                         % hPlot = plot([colV; colV(1)], [rowV; rowV(1)], 'rs-', 'Color', [1 0 0], 'LineStyle', '-');
%                         %                 set(handles.axes5.Children(2), 'XData', colV, 'YData', rowV, 'Color', [1 0 0], 'linewidth', 1,...
%                         %                     'visible', 'on')
%                         hL2 = plot(colV, rowV, 'Color', [0 1 1], 'linewidth', round(get(handles.Slider_zoom, 'Value')));
%                         refreshdata(handles.axes5)
%
%                         %     % GS's structure
%                         %     pointsGS = handles.planC{handles.structfieldnum}(4*(index_selected - 1) + 13).contour(handles.currentSlice_ax).segments(indexseg).points;
%                         %     yCoordsGS = [yCoordsGS; pointsGS(:,2); pointsGS(1,2)]; % + yT; Rotation component (??)  /Users/concetta/OneDrive - Cardiff
%                         %     % University/Concetta_work_cardiff/Projects/Other_projects/CERR-master/CERR_core/Viewers
%                         %     % showStructures line 227
%                         %     xCoordsSG = [xCoordsSG; pointsGS(:,1); pointsGS(1,1)]; %  + xT; Rotation component (??)  /Users/concetta/OneDrive - Cardiff
%                         %     % University/Concetta_work_cardiff/Projects/Other_projects/CERR-master/CERR_core/Viewers
%                         %     % showStructures line 227
%                     end
%                 end
%
% %                 % Legend for the contour
% %                 if isempty(hL2)
% %
% %                 else
% %                     hhl = [hhl hL2];
% %                     labelhhl = [labelhhl {'User volume'}];
% %                 end
%             end
%
% %             if size(hhl, 2) >= 1
% %                 hhllegend = legend(hhl, labelhhl, 'FontSize', 12, 'TextColor', 'white', 'Color', 'Black');
% %             else
% %                 legend(handles.axes5,'hide')
% %             end
%         end
%     case 5
%         if isempty(index_selected) % When no structure is selected
% %             legend(handles.axes5,'hide')
%         else
%             %             % Legend for the contour
%             %             if isempty(hL1)
%             %
%             %             else
%             %             end
%         end
% end

% %% Coronal view
%
% % Setting the axes where visualize the contours
% axes(handles.axes6)
%
% switch n
%
%     case 1
%         if isempty(index_selected) % When no structure is selected
%             legend(handles.axes6,'hide')
%         else
%             % GS's structure
%             if isempty(handles.planC{handles.structfieldnum}(GSsstructure).contour)
%                 xCoords = [];
%                 yCoords = [];
%                 %         set(handles.axes6.Children(1), 'XData', xCoords, 'YData', yCoords, 'Color', [0 1 1], 'linewidth', 1,...
%                 %             'visible', 'off')
%                 hhl = hhl;
%             else
%                 hL1 = [];
%                 numSegs = length(handles.planC{handles.structfieldnum}(GSsstructure).contour(handles.currentSlice_ax).segments);
%
%                 %     thisax = handles.axes6;
%                 %     plot2handle = get(thisax, 'UserData');
%
%                 yCoordsGS = [];
%                 xCoordsSG = [];
%                 for indexseg = 1 : numSegs
%                     pointsM = handles.planC{handles.structfieldnum}(GSsstructure).contour(handles.currentSlice_ax).segments(indexseg).points;
%                     if isempty(pointsM) % Where there are no contours to plot
%                         xCoords = [];
%                         yCoords = [];
%                         %                 set(handles.axes6.Children(1), 'XData', [], 'YData', [], 'Color', [0 1 1], 'linewidth', 1,...
%                         %                     'visible', 'off')
%                         hhl = hhl;
%
%                     else
%                         %Transforming the axial contour to coronal contour
%                         pointsMsag = squeeze(flip(permute(pointsM, [3 2 1]), 1));
%
%                         zValue = handles.planC{handles.scanfieldnum}.scanInfo(handles.currentSlice_ax).zValue;
%
%                         xCoords = pointsMsag(:,1);
%                         yCoords = pointsMsag(:,2);
%
%                         [rowV, colV] = aapmtom(xCoords/scale, yCoords/scale, zValue/scale, xCTOffset/scale, [size(handles.planC{1, 3}.scanArray, 3) imageWidth(2, 1)]);
%
%                         rowV = [rowV; rowV(1)];
%                         colV = [colV; colV(1)];
%
%
%
%                         %Update contour plot
%                         hold on
%                         % hPlot = plot([colV; colV(1)], [rowV; rowV(1)], 'rs-', 'Color', [0 1 1], 'LineStyle', '-');
%                         %                 set(handles.axes6.Children(1), 'XData', colV, 'YData', rowV, 'Color', [0 1 1], 'linewidth', 1,...
%                         %                     'visible', 'on')
%                         hL1 = plot(colV, rowV, 'Color', [1 1 0], 'linewidth', round(get(handles.Slider_zoom, 'Value')));
%                         refreshdata(handles.axes6) %, 'Position', [0.1, 0.012])
%
%                         %                     imshow(blue);
%                         %                 blue = cat(3, zeros(handles.volsize(1 : 2)), zeros(handles.volsize(1 : 2)), ones(handles.volsize(1 : 2)));
%                         %                 set(handles.hShow.Parent.Children(7), 'AlphaData', blue);
%                         %                 https://blogs.mathworks.com/steve/2009/02/18/image-overlay-using-transparency/
%
%                         %     % GS's structure
%                         %     pointsGS = handles.planC{handles.structfieldnum}(4*(index_selected - 1) + 13).contour(handles.currentSlice_ax).segments(indexseg).points;
%                         %     yCoordsGS = [yCoordsGS; pointsGS(:,2); pointsGS(1,2)]; % + yT; Rotation component (??)  /Users/concetta/OneDrive - Cardiff
%                         %     % University/Concetta_work_cardiff/Projects/Other_projects/CERR-master/CERR_core/Viewers
%                         %     % showStructures line 227
%                         %     xCoordsSG = [xCoordsSG; pointsGS(:,1); pointsGS(1,1)]; %  + xT; Rotation component (??)  /Users/concetta/OneDrive - Cardiff
%                         %     % University/Concetta_work_cardiff/Projects/Other_projects/CERR-master/CERR_core/Viewers
%                         %     % showStructures line 227
%                     end
%                 end
%
%                 % Legend for the contour
%                 if isempty(hL1)
%
%                 else
%                     hhl = [hhl hL1];
%                     labelhhl = [labelhhl {'Under region'}];
%                 end
%
%             end
%
%             % User's structure
%             if isempty(handles.planC{handles.structfieldnum}(USsstructure).contour)
%                 xCoords = [];
%                 yCoords = [];
%                 %         set(handles.axes6.Children(2), 'XData', xCoords, 'YData', yCoords, 'Color', [1 0 0], 'linewidth', 1,...
%                 %             'visible', 'off')
%                 hhl = hhl;
%             else
%                 hL2 = [];
%                 numSegs = length(handles.planC{handles.structfieldnum}(USsstructure).contour(handles.currentSlice_ax).segments);
%
%                 yCoordsGS = [];
%                 xCoordsSG = [];
%                 rowV = [];
%                 colV = [];
%                 for indexseg = 1 : numSegs
%                     pointsM = handles.planC{handles.structfieldnum}(USsstructure).contour(handles.currentSlice_ax).segments(indexseg).points;
%
%                     if isempty(pointsM) % Where there are no contours to plot
%                         xCoords = [];
%                         yCoords = [];
%                         %                 set(handles.axes6.Children(2), 'XData', xCoords, 'YData', yCoords, 'Color', [1 0 0], 'linewidth', 1,...
%                         %                     'visible', 'off')
%                         hhl = hhl;
%
%                     else
%                         %Transforming the axial contour to coronal contour
%                         pointsMsag = squeeze(flip(permute(pointsM, [3 2 1]), 1));
%
%                         zValue = handles.planC{handles.scanfieldnum}.scanInfo(handles.currentSlice_ax).zValue;
%
%                         xCoords = pointsMsag(:,1);
%                         yCoords = pointsMsag(:,2);
%
%                         [rowV, colV] = aapmtom(xCoords/scale, yCoords/scale, zValue/scale, xCTOffset/scale, [size(handles.planC{1, 3}.scanArray, 3) imageWidth(2, 1)]);
%
%                         rowV = [rowV; rowV(1)];
%                         colV = [colV; colV(1)];
%
%                         %Update contour plot
%                         hold on
%
%                         % hPlot = plot([colV; colV(1)], [rowV; rowV(1)], 'rs-', 'Color', [1 0 0], 'LineStyle', '-');
%                         %                 set(handles.axes6.Children(2), 'XData', colV, 'YData', rowV, 'Color', [1 0 0], 'linewidth', 1,...
%                         %                     'visible', 'on')
%
%                         hL2 = plot(colV, rowV, 'Color', [0 1 1], 'linewidth', round(get(handles.Slider_zoom, 'Value')));
%                         refreshdata(handles.axes6)
%
%                         %     % GS's structure
%                         %     pointsGS = handles.planC{handles.structfieldnum}(4*(index_selected - 1) + 13).contour(handles.currentSlice_ax).segments(indexseg).points;
%                         %     yCoordsGS = [yCoordsGS; pointsGS(:,2); pointsGS(1,2)]; % + yT; Rotation component (??)  /Users/concetta/OneDrive - Cardiff
%                         %     % University/Concetta_work_cardiff/Projects/Other_projects/CERR-master/CERR_core/Viewers
%                         %     % showStructures line 227
%                         %     xCoordsSG = [xCoordsSG; pointsGS(:,1); pointsGS(1,1)]; %  + xT; Rotation component (??)  /Users/concetta/OneDrive - Cardiff
%                         %     % University/Concetta_work_cardiff/Projects/Other_projects/CERR-master/CERR_core/Viewers
%                         %     % showStructures line 227
%                     end
%                 end
%
%                 % Legend for the contour
%                 if isempty(hL2)
%
%                 else
%                     hhl = [hhl hL2];
%                     labelhhl = [labelhhl {'Over region'}];
%                 end
%             end
%
%             % Structures obtained with the intersection of GS and user's structures
%             if isempty(handles.planC{handles.structfieldnum}(Intersectedstructures).contour)
%                 xCoords = [];
%                 yCoords = [];
%                 %         set(handles.axes6.Children(3), 'XData', xCoords, 'YData', yCoords, 'Color', [0 1 0], 'linewidth', 1,...
%                 %             'visible', 'off')
%                 hhl = hhl;
%             else
%                 hL3 = [];
%                 numSegs = length(handles.planC{handles.structfieldnum}(Intersectedstructures).contour(handles.currentSlice_ax).segments);
%
%                 yCoordsGS = [];
%                 xCoordsSG = [];
%                 rowV = [];
%                 colV = [];
%                 for indexseg = 1 : numSegs
%                     pointsM = handles.planC{handles.structfieldnum}(Intersectedstructures).contour(handles.currentSlice_ax).segments(indexseg).points;
%
%                     if isempty(pointsM) % Where there are no contours to plot
%                         xCoords = [];
%                         yCoords = [];
%                         %                 set(handles.axes6.Children(3), 'XData', xCoords, 'YData', yCoords, 'Color', [0 1 0], 'linewidth', 6,...
%                         %                     'visible', 'off')
%                         hhl = hhl;
%
%                     else
%                         %Transforming the axial contour to coronal contour
%                         pointsMsag = squeeze(flip(permute(pointsM, [3 2 1]), 1));
%
%                         zValue = handles.planC{handles.scanfieldnum}.scanInfo(handles.currentSlice_ax).zValue;
%
%                         xCoords = pointsMsag(:,1);
%                         yCoords = pointsMsag(:,2);
%
%                         [rowV, colV] = aapmtom(xCoords/scale, yCoords/scale, zValue/scale, xCTOffset/scale, [size(handles.planC{1, 3}.scanArray, 3) imageWidth(2, 1)]);
%
%                         rowV = [rowV; rowV(1)];
%                         colV = [colV; colV(1)];
%
%                         %Update contour plot
%                         hold on
%                         %             hPlot = plot([colV; colV(1)], [rowV; rowV(1)], 'rs-', 'Color', [0 1 0], 'LineStyle', '-');
%                         %                 set(handles.axes6.Children(3), 'XData', colV, 'YData', rowV, 'Color', [0 1 0], 'linewidth', 1,...
%                         %                     'visible', 'on')
%                         hL3 = plot(colV, rowV, 'Color', [0 1 0], 'linewidth', round(get(handles.Slider_zoom, 'Value')));
%
%                         %     % GS's structure
%                         %     pointsGS = handles.planC{handles.structfieldnum}(4*(index_selected - 1) + 13).contour(handles.currentSlice_ax).segments(indexseg).points;
%                         %     yCoordsGS = [yCoordsGS; pointsGS(:,2); pointsGS(1,2)]; % + yT; Rotation component (??)  /Users/concetta/OneDrive - Cardiff
%                         %     % University/Concetta_work_cardiff/Projects/Other_projects/CERR-master/CERR_core/Viewers
%                         %     % showStructures line 227
%                         %     xCoordsSG = [xCoordsSG; pointsGS(:,1); pointsGS(1,1)]; %  + xT; Rotation component (??)  /Users/concetta/OneDrive - Cardiff
%                         %     % University/Concetta_work_cardiff/Projects/Other_projects/CERR-master/CERR_core/Viewers
%                         %     % showStructures line 227
%                         refreshdata(handles.axes6)
%                         hold on
%                     end
%                 end
%
%                 % Legend for the contour
%                 if isempty(hL3)
%
%                 else
%                     hhl = [hhl hL3];
%                     labelhhl = [labelhhl {'Common region'}];
%                 end
%             end
%
%             %             if size(hhl, 2) == 3
%             %                 hh = legend(hhl, {'Under region', 'Over region', 'Common region'}, 'FontSize', 12, 'TextColor', 'white', 'Color', 'Black');
%             %             end
%             %             if size(hhl, 2) == 2
%             %                 hh = legend(hhl, {'Under region', 'Over region'}, 'FontSize', 12, 'TextColor', 'white', 'Color', 'Black');
%             %             end
%             %             if size(hhl, 2) == 1
%             %                 if labelhhl(1, 1) == 1
%             %                     hh = legend(hhl, {'Under region'}, 'FontSize', 12, 'TextColor', 'white', 'Color', 'Black');
%             %                 else
%             %                     hh = legend(hhl, {'Over region'}, 'FontSize', 12, 'TextColor', 'white', 'Color', 'Black');
%             %                 end
%             %             end
%             %             if size(hhl, 2) < 1
%             %                 legend(handles.axes6,'hide')
%             %             end
%             if size(hhl, 2) >= 1
%                 hhllegend = legend(hhl, labelhhl, 'FontSize', 12, 'TextColor', 'white', 'Color', 'Black');
%             else
%                 legend(handles.axes6,'hide')
%             end
%         end
%
%     case 2
%         if isempty(index_selected) % When no structure is selected
%             legend(handles.axes6,'hide')
%         else
%             struct_name_case = {'GTV', 'CTVB'}; %, 'PTV'};
%             % User's structure
%             struct_name_temp = struct_name_case{index_selected};
%             % Retrieving information about the selected structure
%             USsstructure  = getStructNum(struct_name_temp, handles.planC, indexS);
%
%
%             if isempty(handles.planC{handles.structfieldnum}(USsstructure).contour)
%                 xCoords = [];
%                 yCoords = [];
%                 %         set(handles.axes6.Children(2), 'XData', xCoords, 'YData', yCoords, 'Color', [1 0 0], 'linewidth', 1,...
%                 %             'visible', 'off')
%                 hhl = hhl;
%             else
%                 hL2 = [];
%                 numSegs = length(handles.planC{handles.structfieldnum}(USsstructure).contour(handles.currentSlice_ax).segments);
%
%                 yCoordsGS = [];
%                 xCoordsSG = [];
%                 rowV = [];
%                 colV = [];
%                 for indexseg = 1 : numSegs
%                     pointsM = handles.planC{handles.structfieldnum}(USsstructure).contour(handles.currentSlice_ax).segments(indexseg).points;
%
%                     if isempty(pointsM) % Where there are no contours to plot
%                         xCoords = [];
%                         yCoords = [];
%                         %                 set(handles.axes6.Children(2), 'XData', xCoords, 'YData', yCoords, 'Color', [1 0 0], 'linewidth', 1,...
%                         %                     'visible', 'off')
%                         hhl = hhl;
%
%                     else
%                         %Transforming the axial contour to coronal contour
%                         pointsMsag = squeeze(flip(permute(pointsM, [3 2 1]), 1));
%
%                         zValue = handles.planC{handles.scanfieldnum}.scanInfo(handles.currentSlice_ax).zValue;
%
%                         xCoords = pointsMsag(:,1);
%                         yCoords = pointsMsag(:,2);
%
%                         [rowV, colV] = aapmtom(xCoords/scale, yCoords/scale, zValue/scale, xCTOffset/scale, [size(handles.planC{1, 3}.scanArray, 3) imageWidth(2, 1)]);
%
%                         rowV = [rowV; rowV(1)];
%                         colV = [colV; colV(1)];
%
%                         %Update contour plot
%                         hold on
%
%                         % hPlot = plot([colV; colV(1)], [rowV; rowV(1)], 'rs-', 'Color', [1 0 0], 'LineStyle', '-');
%                         %                 set(handles.axes6.Children(2), 'XData', colV, 'YData', rowV, 'Color', [1 0 0], 'linewidth', 1,...
%                         %                     'visible', 'on')
%                         hL2 = plot(colV, rowV, 'Color', [0 1 1], 'linewidth', round(get(handles.Slider_zoom, 'Value')));
%                         refreshdata(handles.axes6)
%
%
%                         %     % GS's structure
%                         %     pointsGS = handles.planC{handles.structfieldnum}(4*(index_selected - 1) + 13).contour(handles.currentSlice_ax).segments(indexseg).points;
%                         %     yCoordsGS = [yCoordsGS; pointsGS(:,2); pointsGS(1,2)]; % + yT; Rotation component (??)  /Users/concetta/OneDrive - Cardiff
%                         %     % University/Concetta_work_cardiff/Projects/Other_projects/CERR-master/CERR_core/Viewers
%                         %     % showStructures line 227
%                         %     xCoordsSG = [xCoordsSG; pointsGS(:,1); pointsGS(1,1)]; %  + xT; Rotation component (??)  /Users/concetta/OneDrive - Cardiff
%                         %     % University/Concetta_work_cardiff/Projects/Other_projects/CERR-master/CERR_core/Viewers
%                         %     % showStructures line 227
%                     end
%                 end
%
%                 if isempty(hL2)
%
%                 else
%                     hhl = [hhl, hL2];
%                     labelhhl = [labelhhl, {'User volume'}];
%                 end
%             end
%
%
%         end
%
%         Case = 1; % To be changed based on the analysed case
%         indexS = handles.planC{end};
%         % To be changed based on the analysed case
%         if Case == 1
%             struct_name_case_OARs = {'Vertebra_GS'; 'Aorta_GS'; 'Right lung_GS'; 'Pericardium/great vessels_GS'; ...
%                 'Liver_GS'; 'Stomach_GS'; 'Azygous vein_GS'; 'Left main bronchus_GS'; 'Left lung_GS'};
%
%         end
%         for indexOARs = 1 : size(struct_name_case_OARs, 1) % Retrieving information about the selected structure
%             struct_name_temp = [struct_name_case_OARs{indexOARs}];
%             OARsstructure  = getStructNum(struct_name_temp, handles.planC, indexS);
%
%             % Structures obtained with the intersection of GS and user's structures
%             if isempty(handles.planC{handles.structfieldnum}(OARsstructure).contour)
%                 xCoords = [];
%                 yCoords = [];
%                 %         set(handles.axes6.Children(3), 'XData', xCoords, 'YData', yCoords, 'Color', [0 1 0], 'linewidth', 1,...
%                 %             'visible', 'off')
%                 hhl = hhl;
%             else
%                 hL4 = [];
%                 numSegs = length(handles.planC{handles.structfieldnum}(OARsstructure).contour(handles.currentSlice_ax).segments);
%
%                 yCoordsGS = [];
%                 xCoordsSG = [];
%                 rowV = [];
%                 colV = [];
%                 for indexseg = 1 : numSegs
%                     pointsM = handles.planC{handles.structfieldnum}(OARsstructure).contour(handles.currentSlice_ax).segments(indexseg).points;
%
%                     if isempty(pointsM) % Where there are no contours to plot
%                         xCoords = [];
%                         yCoords = [];
%                         %                 set(handles.axes6.Children(3), 'XData', xCoords, 'YData', yCoords, 'Color', [0 1 0], 'linewidth', 6,...
%                         %                     'visible', 'off')
%                         hhl = hhl;
%
%                     else
%                         %Transforming the axial contour to coronal contour
%                         pointsMsag = squeeze(flip(permute(pointsM, [3 2 1]), 1));
%
%                         zValue = handles.planC{handles.scanfieldnum}.scanInfo(handles.currentSlice_ax).zValue;
%
%                         xCoords = pointsMsag(:,1);
%                         yCoords = pointsMsag(:,2);
%
%                         [rowV, colV] = aapmtom(xCoords/scale, yCoords/scale, zValue/scale, xCTOffset/scale, [size(handles.planC{1, 3}.scanArray, 3) imageWidth(2, 1)]);
%
%                         rowV = [rowV; rowV(1)];
%                         colV = [colV; colV(1)];
%
%                         %Update contour plot
%                         hold on
%                         %             hPlot = plot([colV; colV(1)], [rowV; rowV(1)], 'rs-', 'Color', [0 1 0], 'LineStyle', '-');
%                         %                 set(handles.axes6.Children(3), 'XData', colV, 'YData', rowV, 'Color', [0 1 0], 'linewidth', 1,...
%                         %                     'visible', 'on')
%
%
%
%                         hL4 = plot(colV, rowV, 'color', handles.colorlab(indexOARs, :), 'linewidth', round(get(handles.Slider_zoom, 'Value')));
%
%                         %     % GS's structure
%                         %     pointsGS = handles.planC{handles.structfieldnum}(4*(index_selected - 1) + 13).contour(handles.currentSlice_ax).segments(indexseg).points;
%                         %     yCoordsGS = [yCoordsGS; pointsGS(:,2); pointsGS(1,2)]; % + yT; Rotation component (??)  /Users/concetta/OneDrive - Cardiff
%                         %     % University/Concetta_work_cardiff/Projects/Other_projects/CERR-master/CERR_core/Viewers
%                         %     % showStructures line 227
%                         %     xCoordsSG = [xCoordsSG; pointsGS(:,1); pointsGS(1,1)]; %  + xT; Rotation component (??)  /Users/concetta/OneDrive - Cardiff
%                         %     % University/Concetta_work_cardiff/Projects/Other_projects/CERR-master/CERR_core/Viewers
%                         %     % showStructures line 227
%                         refreshdata(handles.axes6)
%                         hold on
%                     end
%                 end
%
%                 % Legend for the contour
%                 % Only for the first segment. We don't need a label for each segment
%                 if isempty(hL4)
%
%                 else
%                     hhl = [hhl, hL4];
%                     struct_name_temp = struct_name_case_OARs{indexOARs};
%                     labelhhl = [labelhhl {struct_name_temp(1 : end - 3)}];
%                 end
%             end
%         end
%
%         % To add two blank lines in the legend
%         for indexx1 = 1 : 2
%             hL5 = plot([0.1 0.1], [0.1 0.1],'k','LineWidth', 0.1);
%             hhl = [hhl, hL5];
%             struct_name_temp =  ' ';
%             labelhhl = [labelhhl struct_name_temp];
%         end
%
%         % Checking if there is an overlapping between the user selected
%         % structure and each OAR in the current slice. If yes, a flag appears.
%         clear CC
%         S = [];
%         %         new_handle = copyobj(hhllegend, handles.output);
%         %         delete(get(new_handle,'Children'))
%         if isempty(index_selected)
%
%         else
% %             for indexoarmask = 1 : size(handles.OARs_finalmask, 2)
% %                 tempmask = handles.OARs_finalmask{index_selected, indexoarmask};
% %                 flag_label = 0; % flag to avoid adding in the legend two centroids derived from the same OAR
% %                 % Checking if there is an overlapping
% %                 CC(index_selected, indexoarmask) = bwconncomp(tempmask(:, :, handles.currentSlice_ax)); % connectivity = 8 (default)
% %                 if CC(index_selected, indexoarmask).NumObjects > 0
% %                     TT = regionprops(CC(index_selected, indexoarmask), 'Centroid', 'Area');
% %                     S = [S; TT];
% %                     for indexNumObjects = 1 : CC(index_selected, indexoarmask).NumObjects
% %                         if TT(indexNumObjects, 1).Area > 30 % We don't want to annotate just a single point
% %                             % hAnnotAxes = findall(gcf,'type','annotation');
% %                             % arrow_coordinates_x = [(S(1, 1).Centroid(1, 1) - 60), TT(indexNumObjects, 1).Centroid(1, 1)];
% %                             % arrow_coordinates_y = [(S(1, 1).Centroid(1, 2) - 60), TT(indexNumObjects, 1).Centroid(1, 2)];
% %                             % Plots the centroid of the intersection between the US selected structure and the hidden OARs
% %                             hLmask_new = plot(TT(indexNumObjects, 1).Centroid(1, 1), TT(indexNumObjects, 1).Centroid(1, 2), 'Marker', '.', 'color', handles.colorlab(indexoarmask, :),'MarkerSize', 20);
% %                             % from data space to normalized figure coordinates
% %                             % [XFIG, YFIG, DEEP] = ds2fig(arrow_coordinates_x, arrow_coordinates_y);
% %                             % if ~isempty(hAnnotAxes)
% %                             %h_annotation = annotation(handles.output, 'textarrow', XFIG, YFIG);
% %                             % else
% %                             % h_annotation = annotation(handles.output, 'textarrow', XFIG, YFIG, 'String', ['You have included .... in your GTV/CTVB/PTV']);
% %                             % end
% %                             % h_annotation.Color = 'red';
% %                             % h_annotation.FontSize = 14;
% %
% %                             if flag_label == 0 % so to avoid adding in the legend two centroids derived from the same OAR
% %                                 hhl = [hhl, hLmask_new];
% %                                 struct_name_temp = [struct_name_case_OARs{indexoarmask}];
% %                                 labelhhl = [labelhhl , {['You have included too much ' struct_name_temp(1 : end - 3) ' in your ' struct_name_case{index_selected}]}];
% %                                 flag_label = 1;
% %                             end
% %                         else
% %
% %                         end
% %                         % text(round(S(indexNumObjects, 1).Centroid(1, 1)), round(S(indexNumObjects, 1).Centroid(1, 2)),'\rightarrow ciao polla')
% %                     end
% %                 else
% %                     % Nothing appears
% %                 end
% %             end
%
%
% %
% %             if size(hhl, 2) >= 1
% %                 hhllegend = legend(hhl, labelhhl, 'FontSize', 12, 'TextColor', 'white', 'Color', 'Black');
% %             else
% %                 legend(handles.axes6,'hide')
% %             end
%
%             % a2 = axes;
%
%             % %             if size(hhlmask, 2) >= 1
%             % %                 a = axes('position', get(gca,'position'), 'visible', 'on');
%             % %                 hhllegendmask = legend(a, hLmask_new, labelhhlmask, 'FontSize', 12, 'TextColor', 'white', 'Color', 'Black', 'Location', 'NorthWest');
%             % %                 legend update
%             % %             else
%             % %                 a = axes('position',get(gca,'position'),'visible','on');
%             % %                 % legend(a,'hide')
%             % %                 legend update
%             % %             end
%         end
%
%     case 3
%         if isempty(index_selected) % When no structure is selected
%             legend(handles.axes6,'hide')
%         else
%             % GS's structure
%             if isempty(handles.planC{handles.structfieldnum}(GSsstructure_original).contour)
%                 xCoords = [];
%                 yCoords = [];
%                 %         set(handles.axes6.Children(1), 'XData', xCoords, 'YData', yCoords, 'Color', [0 1 1], 'linewidth', 1,...
%                 %             'visible', 'off')
%                 hhl = hhl;
%             else
%                 hL1 = [];
%                 numSegs = length(handles.planC{handles.structfieldnum}(GSsstructure_original).contour(handles.currentSlice_ax).segments);
%
%                 %     thisax = handles.axes6;
%                 %     plot2handle = get(thisax, 'UserData');
%
%                 yCoordsGS = [];
%                 xCoordsSG = [];
%                 for indexseg = 1 : numSegs
%                     pointsM = handles.planC{handles.structfieldnum}(GSsstructure_original).contour(handles.currentSlice_ax).segments(indexseg).points;
%                     if isempty(pointsM) % Where there are no contours to plot
%                         xCoords = [];
%                         yCoords = [];
%                         %                 set(handles.axes6.Children(1), 'XData', [], 'YData', [], 'Color', [0 1 1], 'linewidth', 1,...
%                         %                     'visible', 'off')
%                         hhl = hhl;
%
%                     else
%                         %Transforming the axial contour to coronal contour
%                         pointsMsag = squeeze(flip(permute(pointsM, [3 2 1]), 1));
%
%                         zValue = handles.planC{handles.scanfieldnum}.scanInfo(handles.currentSlice_ax).zValue;
%
%                         xCoords = pointsMsag(:,1);
%                         yCoords = pointsMsag(:,2);
%
%                         [rowV, colV] = aapmtom(xCoords/scale, yCoords/scale, zValue/scale, xCTOffset/scale, [size(handles.planC{1, 3}.scanArray, 3) imageWidth(2, 1)]);
%
%                         rowV = [rowV; rowV(1)];
%                         colV = [colV; colV(1)];
%
%                         %Update contour plot
%                         hold on
%                         % hPlot = plot([colV; colV(1)], [rowV; rowV(1)], 'rs-', 'Color', [0 1 1], 'LineStyle', '-');
%                         %                 set(handles.axes6.Children(1), 'XData', colV, 'YData', rowV, 'Color', [0 1 1], 'linewidth', 1,...
%                         %                     'visible', 'on')
%                         hold on
%                         hL1 = plot(colV, rowV, 'Color', [1 1 0], 'linewidth', round(get(handles.Slider_zoom, 'Value')));
%
%                         refreshdata(handles.axes6) %, 'Position', [0.1, 0.012])
%
%                         %                     imshow(blue);
%                         %                 blue = cat(3, zeros(handles.volsize(1 : 2)), zeros(handles.volsize(1 : 2)), ones(handles.volsize(1 : 2)));
%                         %                 set(handles.hShow.Parent.Children(7), 'AlphaData', blue);
%                         %                 https://blogs.mathworks.com/steve/2009/02/18/image-overlay-using-transparency/
%
%                         %     % GS's structure
%                         %     pointsGS = handles.planC{handles.structfieldnum}(4*(index_selected - 1) + 13).contour(handles.currentSlice_ax).segments(indexseg).points;
%                         %     yCoordsGS = [yCoordsGS; pointsGS(:,2); pointsGS(1,2)]; % + yT; Rotation component (??)  /Users/concetta/OneDrive - Cardiff
%                         %     % University/Concetta_work_cardiff/Projects/Other_projects/CERR-master/CERR_core/Viewers
%                         %     % showStructures line 227
%                         %     xCoordsSG = [xCoordsSG; pointsGS(:,1); pointsGS(1,1)]; %  + xT; Rotation component (??)  /Users/concetta/OneDrive - Cardiff
%                         %     % University/Concetta_work_cardiff/Projects/Other_projects/CERR-master/CERR_core/Viewers
%                         %     % showStructures line 227
%                     end
%                 end
%
%                 % Legend for the contour
%                 if isempty(hL1)
%
%                 else
%                     hhl = [hhl hL1];
%                     labelhhl = [labelhhl {'Reference volume'}];
%                 end
%
%             end
%
%             % User's structure
%             struct_name_temp = struct_name_case{index_selected};
%             % Retrieving information about the selected structure
%             USsstructure  = getStructNum(struct_name_temp, handles.planC, indexS);
%
%
%             if isempty(handles.planC{handles.structfieldnum}(USsstructure).contour)
%                 xCoords = [];
%                 yCoords = [];
%                 %         set(handles.axes6.Children(2), 'XData', xCoords, 'YData', yCoords, 'Color', [1 0 0], 'linewidth', 1,...
%                 %             'visible', 'off')
%                 hhl = hhl;
%             else
%                 hL2 = [];
%                 numSegs = length(handles.planC{handles.structfieldnum}(USsstructure).contour(handles.currentSlice_ax).segments);
%
%                 yCoordsGS = [];
%                 xCoordsSG = [];
%                 rowV = [];
%                 colV = [];
%                 for indexseg = 1 : numSegs
%                     pointsM = handles.planC{handles.structfieldnum}(USsstructure).contour(handles.currentSlice_ax).segments(indexseg).points;
%
%                     if isempty(pointsM) % Where there are no contours to plot
%                         xCoords = [];
%                         yCoords = [];
%                         %                 set(handles.axes6.Children(2), 'XData', xCoords, 'YData', yCoords, 'Color', [1 0 0], 'linewidth', 1,...
%                         %                     'visible', 'off')
%                         hhl = hhl;
%
%                     else
%                         %Transforming the axial contour to coronal contour
%                         pointsMsag = squeeze(flip(permute(pointsM, [3 2 1]), 1));
%
%                         zValue = handles.planC{handles.scanfieldnum}.scanInfo(handles.currentSlice_ax).zValue;
%
%                         xCoords = pointsMsag(:,1);
%                         yCoords = pointsMsag(:,2);
%
%                         [rowV, colV] = aapmtom(xCoords/scale, yCoords/scale, zValue/scale, xCTOffset/scale, [size(handles.planC{1, 3}.scanArray, 3) imageWidth(2, 1)]);
%
%                         rowV = [rowV; rowV(1)];
%                         colV = [colV; colV(1)];
%
%                         %Update contour plot
%                         hold on
%
%                         % hPlot = plot([colV; colV(1)], [rowV; rowV(1)], 'rs-', 'Color', [1 0 0], 'LineStyle', '-');
%                         %                 set(handles.axes6.Children(2), 'XData', colV, 'YData', rowV, 'Color', [1 0 0], 'linewidth', 1,...
%                         %                     'visible', 'on')
%                         hL2 = plot(colV, rowV, 'Color', [0 1 1], 'linewidth', round(get(handles.Slider_zoom, 'Value')));
%                         refreshdata(handles.axes6)
%
%                         %     % GS's structure
%                         %     pointsGS = handles.planC{handles.structfieldnum}(4*(index_selected - 1) + 13).contour(handles.currentSlice_ax).segments(indexseg).points;
%                         %     yCoordsGS = [yCoordsGS; pointsGS(:,2); pointsGS(1,2)]; % + yT; Rotation component (??)  /Users/concetta/OneDrive - Cardiff
%                         %     % University/Concetta_work_cardiff/Projects/Other_projects/CERR-master/CERR_core/Viewers
%                         %     % showStructures line 227
%                         %     xCoordsSG = [xCoordsSG; pointsGS(:,1); pointsGS(1,1)]; %  + xT; Rotation component (??)  /Users/concetta/OneDrive - Cardiff
%                         %     % University/Concetta_work_cardiff/Projects/Other_projects/CERR-master/CERR_core/Viewers
%                         %     % showStructures line 227
%                     end
%                 end
%
%                 % Legend for the contour
%                 if isempty(hL2)
%
%                 else
%                     hhl = [hhl hL2];
%                     labelhhl = [labelhhl {'User volume'}];
%                 end
%             end
%
%             %             % Structures obtained with the intersection of GS and user's structures
%             %             if isempty(handles.planC{handles.structfieldnum}(Intersectedstructures).contour)
%             %                 xCoords = [];
%             %                 yCoords = [];
%             %                 %         set(handles.axes6.Children(3), 'XData', xCoords, 'YData', yCoords, 'Color', [0 1 0], 'linewidth', 1,...
%             %                 %             'visible', 'off')
%             %                 hhl = hhl;
%             %             else
%             %                 hL3 = [];
%             %                 numSegs = length(handles.planC{handles.structfieldnum}(Intersectedstructures).contour(handles.currentSlice_ax).segments);
%             %
%             %                 yCoordsGS = [];
%             %                 xCoordsSG = [];
%             %                 rowV = [];
%             %                 colV = [];
%             %                 for indexseg = 1 : numSegs
%             %                     pointsM = handles.planC{handles.structfieldnum}(Intersectedstructures).contour(handles.currentSlice_ax).segments(indexseg).points;
%             %
%             %                     if isempty(pointsM) % Where there are no contours to plot
%             %                         xCoords = [];
%             %                         yCoords = [];
%             %                         %                 set(handles.axes6.Children(3), 'XData', xCoords, 'YData', yCoords, 'Color', [0 1 0], 'linewidth', 6,...
%             %                         %                     'visible', 'off')
%             %                         hhl = hhl;
%             %
%             %                     else
%             %                         xCoords = pointsM(:,1);
%             %                         yCoords = pointsM(:,2);
%             %
%             %                         [rowV, colV] = aapmtom(xCoords/scale, yCoords/scale, xCTOffset/scale, yCTOffset/scale, imageWidth);
%             %
%             %                         rowV = [rowV; rowV(1)];
%             %                         colV = [colV; colV(1)];
%             %
%             %                         %Update contour plot
%             %                         hold on
%             %                         %             hPlot = plot([colV; colV(1)], [rowV; rowV(1)], 'rs-', 'Color', [0 1 0], 'LineStyle', '-');
%             %                         %                 set(handles.axes6.Children(3), 'XData', colV, 'YData', rowV, 'Color', [0 1 0], 'linewidth', 1,...
%             %                         %                     'visible', 'on')
%             %                         hL3 = plot(colV, rowV, 'Color', [0 1 0], 'linewidth', 1);
%             %
%             %                         %     % GS's structure
%             %                         %     pointsGS = handles.planC{handles.structfieldnum}(4*(index_selected - 1) + 13).contour(handles.currentSlice_ax).segments(indexseg).points;
%             %                         %     yCoordsGS = [yCoordsGS; pointsGS(:,2); pointsGS(1,2)]; % + yT; Rotation component (??)  /Users/concetta/OneDrive - Cardiff
%             %                         %     % University/Concetta_work_cardiff/Projects/Other_projects/CERR-master/CERR_core/Viewers
%             %                         %     % showStructures line 227
%             %                         %     xCoordsSG = [xCoordsSG; pointsGS(:,1); pointsGS(1,1)]; %  + xT; Rotation component (??)  /Users/concetta/OneDrive - Cardiff
%             %                         %     % University/Concetta_work_cardiff/Projects/Other_projects/CERR-master/CERR_core/Viewers
%             %                         %     % showStructures line 227
%             %                         refreshdata(handles.axes6)
%             %                         hold on
%             %                     end
%             %                 end
%             %
%             %                 % Legend for the contour
%             %                 if isempty(hL3)
%             %
%             %                 else
%             %                     hhl = [hhl hL3];
%             %                     labelhhl = [labelhhl {'Common region'}];
%             %                 end
%             %             end
%
%             % GS's structure max
%             if isempty(handles.planC{handles.structfieldnum}(GSsstructuremax).contour)
%                 xCoords = [];
%                 yCoords = [];
%                 %         set(handles.axes6.Children(3), 'XData', xCoords, 'YData', yCoords, 'Color', [0 1 0], 'linewidth', 1,...
%                 %             'visible', 'off')
%                 hhl = hhl;
%             else
%                 hL4 = [];
%                 numSegs = length(handles.planC{handles.structfieldnum}(GSsstructuremax).contour(handles.currentSlice_ax).segments);
%
%                 yCoordsGS = [];
%                 xCoordsSG = [];
%                 rowV = [];
%                 colV = [];
%                 for indexseg = 1 : numSegs
%                     pointsM = handles.planC{handles.structfieldnum}(GSsstructuremax).contour(handles.currentSlice_ax).segments(indexseg).points;
%
%                     if isempty(pointsM) % Where there are no contours to plot
%                         xCoords = [];
%                         yCoords = [];
%                         %                 set(handles.axes6.Children(3), 'XData', xCoords, 'YData', yCoords, 'Color', [0 1 0], 'linewidth', 6,...
%                         %                     'visible', 'off')
%                         hhl = hhl;
%
%                     else
%                         %Transforming the axial contour to coronal contour
%                         pointsMsag = squeeze(flip(permute(pointsM, [3 2 1]), 1));
%
%                         zValue = handles.planC{handles.scanfieldnum}.scanInfo(handles.currentSlice_ax).zValue;
%
%                         xCoords = pointsMsag(:,1);
%                         yCoords = pointsMsag(:,2);
%
%                         [rowV, colV] = aapmtom(xCoords/scale, yCoords/scale, zValue/scale, xCTOffset/scale, [size(handles.planC{1, 3}.scanArray, 3) imageWidth(2, 1)]);
%
%                         rowV = [rowV; rowV(1)];
%                         colV = [colV; colV(1)];
%
%                         %Update contour plot
%                         hold on
%                         %             hPlot = plot([colV; colV(1)], [rowV; rowV(1)], 'rs-', 'Color', [0 1 0], 'LineStyle', '-');
%                         %                 set(handles.axes6.Children(3), 'XData', colV, 'YData', rowV, 'Color', [0 1 0], 'linewidth', 1,...
%                         %                     'visible', 'on')
%                         % hL4 = plot(colV, rowV, 'Color', [0 1 1], 'linewidth', 1, 'LineStyle', '--');
%                         hL4 = plot(colV, rowV, 'Color', [1 0 0], 'linewidth', round(get(handles.Slider_zoom, 'Value')));
%
%                         %     % GS's structure
%                         %     pointsGS = handles.planC{handles.structfieldnum}(4*(index_selected - 1) + 13).contour(handles.currentSlice_ax).segments(indexseg).points;
%                         %     yCoordsGS = [yCoordsGS; pointsGS(:,2); pointsGS(1,2)]; % + yT; Rotation component (??)  /Users/concetta/OneDrive - Cardiff
%                         %     % University/Concetta_work_cardiff/Projects/Other_projects/CERR-master/CERR_core/Viewers
%                         %     % showStructures line 227
%                         %     xCoordsSG = [xCoordsSG; pointsGS(:,1); pointsGS(1,1)]; %  + xT; Rotation component (??)  /Users/concetta/OneDrive - Cardiff
%                         %     % University/Concetta_work_cardiff/Projects/Other_projects/CERR-master/CERR_core/Viewers
%                         %     % showStructures line 227
%                         refreshdata(handles.axes6)
%                         hold on
%                     end
%                 end
%
%                 % Legend for the contour
%                 if isempty(hL4)
%
%                 else
%                     hhl = [hhl hL4];
%                     struct_name_temp = struct_name_case{index_selected};
%                     labelhhl = [labelhhl {[struct_name_temp ' max']}];
%                 end
%             end
%
%             % GS's structure min
%             if isempty(handles.planC{handles.structfieldnum}(GSsstructuremin).contour)
%                 xCoords = [];
%                 yCoords = [];
%                 %         set(handles.axes6.Children(3), 'XData', xCoords, 'YData', yCoords, 'Color', [0 1 0], 'linewidth', 1,...
%                 %             'visible', 'off')
%                 hhl = hhl;
%             else
%                 hL5 = [];
%                 numSegs = length(handles.planC{handles.structfieldnum}(GSsstructuremin).contour(handles.currentSlice_ax).segments);
%
%                 yCoordsGS = [];
%                 xCoordsSG = [];
%                 rowV = [];
%                 colV = [];
%                 for indexseg = 1 : numSegs
%                     pointsM = handles.planC{handles.structfieldnum}(GSsstructuremin).contour(handles.currentSlice_ax).segments(indexseg).points;
%
%                     if isempty(pointsM) % Where there are no contours to plot
%                         xCoords = [];
%                         yCoords = [];
%                         %                 set(handles.axes6.Children(3), 'XData', xCoords, 'YData', yCoords, 'Color', [0 1 0], 'linewidth', 6,...
%                         %                     'visible', 'off')
%                         hhl = hhl;
%
%                     else
%                         %Transforming the axial contour to coronal contour
%                         pointsMsag = squeeze(flip(permute(pointsM, [3 2 1]), 1));
%
%                         zValue = handles.planC{handles.scanfieldnum}.scanInfo(handles.currentSlice_ax).zValue;
%
%                         xCoords = pointsMsag(:,1);
%                         yCoords = pointsMsag(:,2);
%
%                         [rowV, colV] = aapmtom(xCoords/scale, yCoords/scale, zValue/scale, xCTOffset/scale, [size(handles.planC{1, 3}.scanArray, 3) imageWidth(2, 1)]);
%
%                         rowV = [rowV; rowV(1)];
%                         colV = [colV; colV(1)];
%
%                         %Update contour plot
%                         hold on
%                         %             hPlot = plot([colV; colV(1)], [rowV; rowV(1)], 'rs-', 'Color', [0 1 0], 'LineStyle', '-');
%                         %                 set(handles.axes6.Children(3), 'XData', colV, 'YData', rowV, 'Color', [0 1 0], 'linewidth', 1,...
%                         %                     'visible', 'on')
%                         % hL5 = plot(colV, rowV, 'Color', [1 1 0], 'linewidth', 2, 'LineStyle', ':');
%                         hL5 = plot(colV, rowV, 'Color', [0 1 0], 'linewidth', round(get(handles.Slider_zoom, 'Value')));
%
%                         %     % GS's structure
%                         %     pointsGS = handles.planC{handles.structfieldnum}(4*(index_selected - 1) + 13).contour(handles.currentSlice_ax).segments(indexseg).points;
%                         %     yCoordsGS = [yCoordsGS; pointsGS(:,2); pointsGS(1,2)]; % + yT; Rotation component (??)  /Users/concetta/OneDrive - Cardiff
%                         %     % University/Concetta_work_cardiff/Projects/Other_projects/CERR-master/CERR_core/Viewers
%                         %     % showStructures line 227
%                         %     xCoordsSG = [xCoordsSG; pointsGS(:,1); pointsGS(1,1)]; %  + xT; Rotation component (??)  /Users/concetta/OneDrive - Cardiff
%                         %     % University/Concetta_work_cardiff/Projects/Other_projects/CERR-master/CERR_core/Viewers
%                         %     % showStructures line 227
%                         refreshdata(handles.axes6)
%                         hold on
%                     end
%                 end
%
%                 % Legend for the contour
%                 if isempty(hL5)
%
%                 else
%                     hhl = [hhl hL5];
%                     struct_name_temp = struct_name_case{index_selected};
%                     labelhhl = [labelhhl {[struct_name_temp ' min']}];
%                 end
%             end
%
%             %             if size(hhl, 2) == 3
%             %                 hh = legend(hhl, {'Under region', 'Over region', 'Common region'}, 'FontSize', 12, 'TextColor', 'white', 'Color', 'Black');
%             %             end
%             %             if size(hhl, 2) == 2
%             %                 hh = legend(hhl, {'Under region', 'Over region'}, 'FontSize', 12, 'TextColor', 'white', 'Color', 'Black');
%             %             end
%             %             if size(hhl, 2) == 1
%             %                 if labelhhl(1, 1) == 1
%             %                     hh = legend(hhl, {'Under region'}, 'FontSize', 12, 'TextColor', 'white', 'Color', 'Black');
%             %                 else
%             %                     hh = legend(hhl, {'Over region'}, 'FontSize', 12, 'TextColor', 'white', 'Color', 'Black');
%             %                 end
%             %             end
%             %             if size(hhl, 2) < 1
%             %                 legend(handles.axes6,'hide')
%             %             end
%             if size(hhl, 2) >= 1
%                 hhllegend = legend(hhl, labelhhl, 'FontSize', 12, 'TextColor', 'white', 'Color', 'Black');
%             else
%                 legend(handles.axes6,'hide')
%             end
%         end
%
%     case 4
%         if isempty(index_selected) % When no structure is selected
%             legend(handles.axes6,'hide')
%         else
%             % GS's structure
%             if isempty(handles.planC{handles.structfieldnum}(GSsstructure_original).contour)
%                 xCoords = [];
%                 yCoords = [];
%                 %         set(handles.axes6.Children(1), 'XData', xCoords, 'YData', yCoords, 'Color', [0 1 1], 'linewidth', 1,...
%                 %             'visible', 'off')
%                 hhl = hhl;
%             else
%                 hL1 = [];
%                 numSegs = length(handles.planC{handles.structfieldnum}(GSsstructure_original).contour(handles.currentSlice_ax).segments);
%
%                 %     thisax = handles.axes6;
%                 %     plot2handle = get(thisax, 'UserData');
%
%                 yCoordsGS = [];
%                 xCoordsSG = [];
%                 for indexseg = 1 : numSegs
%                     pointsM = handles.planC{handles.structfieldnum}(GSsstructure_original).contour(handles.currentSlice_ax).segments(indexseg).points;
%                     if isempty(pointsM) % Where there are no contours to plot
%                         xCoords = [];
%                         yCoords = [];
%                         %                 set(handles.axes6.Children(1), 'XData', [], 'YData', [], 'Color', [0 1 1], 'linewidth', 1,...
%                         %                     'visible', 'off')
%                         hhl = hhl;
%
%                     else
%                         %Transforming the axial contour to coronal contour
%                         pointsMsag = squeeze(flip(permute(pointsM, [3 2 1]), 1));
%
%                         zValue = handles.planC{handles.scanfieldnum}.scanInfo(handles.currentSlice_ax).zValue;
%
%                         xCoords = pointsMsag(:,1);
%                         yCoords = pointsMsag(:,2);
%
%                         [rowV, colV] = aapmtom(xCoords/scale, yCoords/scale, zValue/scale, xCTOffset/scale, [size(handles.planC{1, 3}.scanArray, 3) imageWidth(2, 1)]);
%
%                         rowV = [rowV; rowV(1)];
%                         colV = [colV; colV(1)];
%
%                         %Update contour plot
%                         % hold on
%                         % hPlot = plot([colV; colV(1)], [rowV; rowV(1)], 'rs-', 'Color', [0 1 1], 'LineStyle', '-');
%                         %                 set(handles.axes6.Children(1), 'XData', colV, 'YData', rowV, 'Color', [0 1 1], 'linewidth', 1,...
%                         %                     'visible', 'on')
%                         hL1 = plot(colV, rowV, 'Color', [0 1 1], 'linewidth', round(get(handles.Slider_zoom, 'Value')));
%                         refreshdata(handles.axes6) %, 'Position', [0.1, 0.012])
%
%                         %                     imshow(blue);
%                         %                 blue = cat(3, zeros(handles.volsize(1 : 2)), zeros(handles.volsize(1 : 2)), ones(handles.volsize(1 : 2)));
%                         %                 set(handles.hShow.Parent.Children(7), 'AlphaData', blue);
%                         %                 https://blogs.mathworks.com/steve/2009/02/18/image-overlay-using-transparency/
%
%                         %     % GS's structure
%                         %     pointsGS = handles.planC{handles.structfieldnum}(4*(index_selected - 1) + 13).contour(handles.currentSlice_ax).segments(indexseg).points;
%                         %     yCoordsGS = [yCoordsGS; pointsGS(:,2); pointsGS(1,2)]; % + yT; Rotation component (??)  /Users/concetta/OneDrive - Cardiff
%                         %     % University/Concetta_work_cardiff/Projects/Other_projects/CERR-master/CERR_core/Viewers
%                         %     % showStructures line 227
%                         %     xCoordsSG = [xCoordsSG; pointsGS(:,1); pointsGS(1,1)]; %  + xT; Rotation component (??)  /Users/concetta/OneDrive - Cardiff
%                         %     % University/Concetta_work_cardiff/Projects/Other_projects/CERR-master/CERR_core/Viewers
%                         %     % showStructures line 227
%                     end
%                 end
%
%                 % Legend for the contour
%                 if isempty(hL1)
%
%                 else
%                     hhl = [hhl hL1];
%                     labelhhl = [labelhhl {'Reference volume'}];
%                 end
%
%             end
%
%             % User's structure
%             if isempty(handles.planC{handles.structfieldnum}(USsstructure_original).contour)
%                 xCoords = [];
%                 yCoords = [];
%                 %         set(handles.axes6.Children(2), 'XData', xCoords, 'YData', yCoords, 'Color', [1 0 0], 'linewidth', 1,...
%                 %             'visible', 'off')
%                 hhl = hhl;
%             else
%                 hL2 = [];
%                 numSegs = length(handles.planC{handles.structfieldnum}(USsstructure_original).contour(handles.currentSlice_ax).segments);
%
%                 yCoordsGS = [];
%                 xCoordsSG = [];
%                 rowV = [];
%                 colV = [];
%                 for indexseg = 1 : numSegs
%                     pointsM = handles.planC{handles.structfieldnum}(USsstructure_original).contour(handles.currentSlice_ax).segments(indexseg).points;
%
%                     if isempty(pointsM) % Where there are no contours to plot
%                         xCoords = [];
%                         yCoords = [];
%                         %                 set(handles.axes6.Children(2), 'XData', xCoords, 'YData', yCoords, 'Color', [1 0 0], 'linewidth', 1,...
%                         %                     'visible', 'off')
%                         hhl = hhl;
%
%                     else
%                         %Transforming the axial contour to coronal contour
%                         pointsMsag = squeeze(flip(permute(pointsM, [3 2 1]), 1));
%
%                         zValue = handles.planC{handles.scanfieldnum}.scanInfo(handles.currentSlice_ax).zValue;
%
%                         xCoords = pointsMsag(:,1);
%                         yCoords = pointsMsag(:,2);
%
%                         [rowV, colV] = aapmtom(xCoords/scale, yCoords/scale, zValue/scale, xCTOffset/scale, [size(handles.planC{1, 3}.scanArray, 3) imageWidth(2, 1)]);
%
%                         rowV = [rowV; rowV(1)];
%                         colV = [colV; colV(1)];
%
%                         %Update contour plot
%                         hold on
%
%                         % hPlot = plot([colV; colV(1)], [rowV; rowV(1)], 'rs-', 'Color', [1 0 0], 'LineStyle', '-');
%                         %                 set(handles.axes6.Children(2), 'XData', colV, 'YData', rowV, 'Color', [1 0 0], 'linewidth', 1,...
%                         %                     'visible', 'on')
%                         hL2 = plot(colV, rowV, 'Color', [0 1 1], 'linewidth', round(get(handles.Slider_zoom, 'Value')));
%                         refreshdata(handles.axes6)
%
%                         %     % GS's structure
%                         %     pointsGS = handles.planC{handles.structfieldnum}(4*(index_selected - 1) + 13).contour(handles.currentSlice_ax).segments(indexseg).points;
%                         %     yCoordsGS = [yCoordsGS; pointsGS(:,2); pointsGS(1,2)]; % + yT; Rotation component (??)  /Users/concetta/OneDrive - Cardiff
%                         %     % University/Concetta_work_cardiff/Projects/Other_projects/CERR-master/CERR_core/Viewers
%                         %     % showStructures line 227
%                         %     xCoordsSG = [xCoordsSG; pointsGS(:,1); pointsGS(1,1)]; %  + xT; Rotation component (??)  /Users/concetta/OneDrive - Cardiff
%                         %     % University/Concetta_work_cardiff/Projects/Other_projects/CERR-master/CERR_core/Viewers
%                         %     % showStructures line 227
%                     end
%                 end
%
%                 % Legend for the contour
%                 if isempty(hL2)
%
%                 else
%                     hhl = [hhl hL2];
%                     labelhhl = [labelhhl {'User volume'}];
%                 end
%             end
%
%             if size(hhl, 2) >= 1
%                 hhllegend = legend(hhl, labelhhl, 'FontSize', 12, 'TextColor', 'white', 'Color', 'Black');
%             else
%                 legend(handles.axes6,'hide')
%             end
%         end
%     case 5
%         if isempty(index_selected) % When no structure is selected
%             legend(handles.axes6,'hide')
%         else
%             %             % Legend for the contour
%             %             if isempty(hL1)
%             %
%             %             else
%             %             end
%         end
% end

end
