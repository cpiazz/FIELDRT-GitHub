function [indeximagestotal, waitingbar_index] = FIELDRTPdffiguresaving(varargin, range_total, struct_name_case, module, struct_name_case_OARs, module_imagetoshow, wb, waitingbar_index, waitingbar_struct_name_case)
% Creating and storing the figure to save in the pdf report

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
% Concetta Piazzese June 2018

% Retrieving the informaion
handles.planC = varargin{1, 1}.planC;
handles.numofStruct = varargin{1, 1}.numofStruct;
handles.FIELDRT = varargin{1, 1}.FIELDRT;
handles.jaccard3Doutput = varargin{1, 1}.jaccard3Doutput;
handles.jaccard2Doutput = varargin{1, 1}.jaccard2Doutput;
handles.volumestat_GS = varargin{1, 1}.volumestat_GS;
handles.volumestat_US = varargin{1, 1}.volumestat_US;
handles.volume_Ratio = varargin{1, 1}.volume_Ratio;
% handles.colorlab = varargin{1, 1}.colorlab;
handles.OARs_finalmask = varargin{1, 1}.OARs_finalmask;
handles.OARsdummy_finalmask = varargin{1, 1}.OARsdummy_finalmask;
handles.Case = varargin{1, 1}.Case;
handles.FIELDRTDataAttempts = varargin{1, 1}.FIELDRTDataAttempts;
handles.AttemptInformation = varargin{1, 1}.AttemptInformation;
handles.FIELDRTGSCases = varargin{1, 1}.FIELDRTGSCases;

% Retrieving some information
indexS = handles.planC{end};
handles.structfieldnum = indexS.structures;
handles.scanfieldnum = indexS.scan;
handles.volsize = size(handles.planC{1, handles.scanfieldnum}.scanArray);

% flag to create the first page of the report with 4 images just one time
% no montage with 4 images needed for the min and max module
if strcmp(module, 'Overunder') || strcmp(module, 'OARs')
    firstpageflag = 1; % flag on --> montage with 4 images will be created
else
    firstpageflag = 0; % flag off --> no montage with 4 images will be created
end

%% Creating the images for each structure and module 
for index1 = 1 : size(struct_name_case, 2)
    % Name of the structure
    struct_name_temp = struct_name_case{index1};
    
    
    % Updating the waiting bar
    wb = waitbar((waitingbar_index / waitingbar_struct_name_case), wb);
    
    % Index for counting the images created in each folder
    % indeximages_total = 0;
    
    % Getting the information about the screen size
    screen = get(0,'ScreenSize');
    
    % figure for the images
    f = figure('visible','off', 'Position', [1915 1915 1110 1150]);
    
    % To have the figure with Minimal White Space
    ax = gca;
    outerpos = ax.OuterPosition;
    ti = ax.TightInset;
    left = outerpos(1) + ti(1);
    bottom = outerpos(2) + ti(2);
    ax_width = outerpos(3) - ti(1) - ti(3);
    ax_height = outerpos(4) - ti(2)*5 - ti(4);
    ax.Position = [left bottom ax_width ax_height];
    
    % To avoid text to becoma white when printing the figure
    set(gcf, 'InvertHardCopy', 'off');
    
    % Retrieving information about the selected structure
    structnum  = getStructNum([struct_name_case{index1} '_GS'], handles.planC, indexS);
    
    
    if strcmp(module, 'Overunder')
        indexmodule = 1;
    end
    
    if strcmp(module, 'OARs')
        indexmodule = 2;
    end
    
    if strcmp(module, 'Minmax')
        indexmodule = 3;
    end
    
    switch module
        case 'Overunder'
            % Index for counting the images created
            indeximages_temp = 0;
            
            for index2 = range_total(index1, 1) : range_total(index1, 2)
                
                if module_imagetoshow(index1, index2) == 1 % The image was flagged so it needs to be shown
                    hhl = [];
                    labelhhl = [];
                    
                    % Range of slices
                    minSlice = min(min(handles.planC{1, handles.scanfieldnum}.scanArray(:, :, index2)));
                    maxSlice = max(max(handles.planC{1, handles.scanfieldnum}.scanArray(:, :, index2)));
                    
                    % Retieving information for contours plotting
                    scale =  handles.planC{handles.scanfieldnum}.scanInfo(index2).grid1Units;
                    imageWidth(1, 1) =  handles.planC{handles.scanfieldnum}.scanInfo(index2).sizeOfDimension1;
                    imageWidth(2, 1) =  handles.planC{handles.scanfieldnum}.scanInfo(index2).sizeOfDimension2;
                    
                    if ~isempty(handles.planC{handles.scanfieldnum}.scanInfo(1).xOffset)
                        xCTOffset = handles.planC{handles.scanfieldnum}.scanInfo(1).xOffset;
                        yCTOffset = handles.planC{handles.scanfieldnum}.scanInfo(1).yOffset;
                    else
                        xCTOffset = 0;
                        yCTOffset = 0;
                    end
                    
                    pointsM = [];
                    
                    % Visualizing the image
                    hhh = imshow (handles.planC{1, handles.scanfieldnum}.scanArray(:, :, index2), [minSlice maxSlice], 'InitialMagnification','fit');
                    
                    hold on
                    
                    %% GS's structure
                    GSsstructure  = getStructNum([struct_name_temp '_GSdiffinters'], handles.planC, indexS);
                    
                    hL1 = [];
                    
                    yCoordsGS = [];
                    xCoordsSG = [];
                    
                    % Checking if there are points in the contour
                    if (isempty(handles.planC{handles.structfieldnum}(GSsstructure).contour(index2).segments))  % Where there are no contours to plot
                        
                        clear hL1
                        
                    else
                        
                        pointsM = handles.planC{handles.structfieldnum}(GSsstructure).contour(index2).segments.points;
                        
                        if isempty(pointsM) % Another check if there are points in the contour
                            
                        else
                            
                            xCoords = pointsM(:,1);
                            yCoords = pointsM(:,2);
                            
                            [rowV, colV] = aapmtom(xCoords/scale, yCoords/scale, xCTOffset/scale, yCTOffset/scale, imageWidth);
                            
                            rowV = [rowV; rowV(1)];
                            colV = [colV; colV(1)];
                            
                            % Plotting the contour
                            hL1 = plot(colV, rowV, 'Color', [0 1 1], 'linewidth', 4);
                        end
                    end
                    
                    % Legend for the contour
                    if exist('hL1')
                        if isempty(hL1)
                            
                        else
                            hhl = [hhl hL1];
                            labelhhl = [labelhhl {'Under region'}];
                        end
                    else
                    end
                    
                    
                    %% User's structure
                    USsstructure  = getStructNum([struct_name_temp '_USdiffinters'], handles.planC, indexS);
                    
                    hL2 = [];
                    
                    yCoordsGS = [];
                    xCoordsSG = [];
                    
                    % Checking if there are points in the contour
                    if isempty(handles.planC{handles.structfieldnum}(USsstructure).contour(index2).segments) % Where there are no contours to plot
                        
                        clear hL2
                        
                    else
                        pointsM = handles.planC{handles.structfieldnum}(USsstructure).contour(index2).segments.points;
                        
                        if isempty(pointsM) % Another check if there are points in the contour
                            
                        else
                            xCoords = pointsM(:,1);
                            yCoords = pointsM(:,2);
                            
                            [rowV, colV] = aapmtom(xCoords/scale, yCoords/scale, xCTOffset/scale, yCTOffset/scale, imageWidth);
                            
                            rowV = [rowV; rowV(1)];
                            colV = [colV; colV(1)];
                            
                            % Plotting the contour
                            hL2 = plot(colV, rowV, 'Color', [1 0 0], 'linewidth', 4);
                        end
                    end
                    
                    % Legend for the contour
                    if exist('hL2')
                        if isempty(hL2)
                            
                        else
                            hhl = [hhl hL2];
                            labelhhl = [labelhhl {'Over region'}];
                        end
                    else
                    end
                    
                    % Structures obtained with the intersection of GS and user's structures
                    Intersectedstructures  = getStructNum([struct_name_temp '_inters'], handles.planC, indexS);
                    
                    hL3 = [];
                    
                    yCoordsGS = [];
                    xCoordsSG = [];
                    
                    % Checking if there are points in the contour
                    if isempty(handles.planC{handles.structfieldnum}(Intersectedstructures).contour(index2).segments) % Where there are no contours to plot
                        
                    else
                        
                        pointsM = handles.planC{handles.structfieldnum}(Intersectedstructures).contour(index2).segments.points;
                        
                        if isempty(pointsM) % Another check if there are points in the contour
                            
                        else
                            xCoords = pointsM(:,1);
                            yCoords = pointsM(:,2);
                            
                            [rowV, colV] = aapmtom(xCoords/scale, yCoords/scale, xCTOffset/scale, yCTOffset/scale, imageWidth);
                            
                            rowV = [rowV; rowV(1)];
                            colV = [colV; colV(1)];
                            
                            % Plotting the contour
                            hL3 = plot(colV, rowV, 'Color', [0 1 0], 'linewidth', 4);
                        end
                    end
                    
                    % Legend for the contour
                    if exist('hL3')
                        if isempty(hL3)
                            
                        else
                            hhl = [hhl hL3];
                            labelhhl = [labelhhl {'Common region region'}];
                        end
                    else
                    end
                    
                    if size(hhl, 2) >= 1
                        hhllegend = legend(hhl, labelhhl, 'FontSize', 32, 'TextColor', 'white', 'Color', 'Black');
                    else
                        
                    end
                    
                    % Title of the image
                    caption = sprintf('Slice %i', index2);
                    title(caption, 'FontSize', 40, 'Color' , 'k');
                    
                    
                    indeximages_temp = indeximages_temp + 1;
                    
                    if ispc
                        % saveas(f, [handles.FIELDRT.env.userDir '\FIELDRT_Reports\Images\' module '\'  struct_name_case{index1} '\' struct_name_temp int2str(index2)] , 'tiff')
                        print(f, [handles.FIELDRT.env.userDir '\FIELDRT_Reports\Images\' module '\'  struct_name_case{index1} '\' struct_name_temp int2str(index2)] , '-dtiffn')
                    else
                        % saveas(f, [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp int2str(index2)] , 'tiff')
                        print(f, [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp int2str(index2)] , '-dtiffn')
                    end
                    
                    hold off
                end
            end
            
            % indeximages_total(1) = indeximages_temp;
            
            % closing the figure
            clear f
            
            
            %% Creating the montage with the created figures
            
            %     % figure for the montage
            %     f = figure('visible','off', 'Position', [screen(1) screen(2) (screen(3) + 500) (screen(4) + 200) ]);
            %
            %     % To have the figure with Minimal White Space
            %     ax = gca;
            %     outerpos = ax.OuterPosition;
            %     ti = ax.TightInset;
            %     left = outerpos(1) + ti(1);
            %     bottom = outerpos(2) + ti(2);
            %     ax_width = outerpos(3) - ti(1) - ti(3);
            %     ax_height = outerpos(4) - ti(2) - ti(4);
            %     ax.Position = [left bottom ax_width ax_height];
            
            
            if indeximages_temp == 0 % no images are created so no montage to create and save
                
                indeximagestotal (1, index1) = 0;
                
            else % there are images so a montage need to be created
                
                nummont = 0;
                
                if firstpageflag == 1 
                    if ispc
                        fileFolder = [handles.FIELDRT.env.userDir '\FIELDRT_Reports\Images\' module '\' struct_name_case{index1} '\'];
                    else
                        fileFolder = [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/'];
                    end
                    
                    dirOutput_iniz = dir(fullfile(fileFolder,[struct_name_case{index1} '*.tif']));
                    
                    % Only for the first structure a collage of four images is
                    % created. This is needed because of the layout of the first page of
                    % each functionality (i.e over/under, OARs...)
                    if size(dirOutput_iniz, 1) < 4
                        dirOutput_temp = dirOutput_iniz(1 : size(dirOutput_iniz, 1));
                        fileNames_temp = string({dirOutput_temp.name}); %retrieving filenames
                        firstpageflag = 0;
                    else
                        dirOutput_temp = dirOutput_iniz(1 : 4);
                        fileNames_temp = string({dirOutput_temp.name}); %retrieving filenames
                        dirOutput = dirOutput_iniz(5 : end);
                        firstpageflag = 0;
                    end
                    
                    % Creation of the first collage with four images for the first
                    % structure
                    clear I
                    
                    if size(dirOutput_temp, 1) > 0
                        for indexfig = 1 : size(dirOutput_temp, 1)
                            if ispc
                                I{indexfig} = imread([dirOutput_temp(indexfig).folder '\' dirOutput_temp(indexfig).name]);
                            else
                                I{indexfig} = imread([dirOutput_temp(indexfig).folder '/' dirOutput_temp(indexfig).name]);
                            end
                        end
                        
                        if size(dirOutput_temp, 1) == 1
                             Im = [I{1, 1} ones(size(I{1, 1}))*255;
                                 ones(size(I{1, 1}))*255 ones(size(I{1, 1}))*255];
                        end
                        
                        if size(dirOutput_temp, 1) == 2
                             Im = [I{1, 1} I{1, 2};
                                 ones(size(I{1, 1}))*255 ones(size(I{1, 1}))*255];
                        end
                        
                        if size(dirOutput_temp, 1) == 3
                             Im = [I{1, 1} I{1, 2};
                            I{1, 3} ones(size(I{1, 1}))*255];
                        end
                        
                        if size(dirOutput_temp, 1) == 4
                             Im = [I{1, 1} I{1, 2};
                            I{1, 3} I{1, 4}];
                        end
                        
                        if ispc
                            % saveas(f, [handles.FIELDRT.env.userDir '\FIELDRT_Reports\Images\' module '\' struct_name_case{index1} '\' struct_name_temp '_montage_' int2str(indexmontage)] , 'tiff')
                            % print(f, [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage)] , '-dtiffn', '-r500')
                            imwrite(Im,[handles.FIELDRT.env.userDir '\FIELDRT_Reports\Images\' module '\' struct_name_case{index1} '\' struct_name_temp '_montage_' int2str(1) '.tiff'],'tiff')
                        else
                            % saveas(f, [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage)] , 'tiff')
                            % print(f, [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage)] , '-dtiffn', '-r500')
                            imwrite(Im,[handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(1) '.tiff'],'tiff')
                        end
                    end
                    
                    
                    % creation of the montage with 6 images
                    if size(dirOutput_iniz,1) > 4 
                        
                        nummont = ceil(size(dirOutput,1) / 6);
                        
                        for indexmontage = 1 : nummont
                            
                            if indexmontage == nummont
                                
                                clear I
                                
                                % number of images not available to fill the 3*2 collage
                                % imadiff = (nummont * 6) - size(dirOutput,1);
                                
                                for indexfig = (6 * (indexmontage - 1) + 1) : size(dirOutput, 1)
                                    if ispc
                                        I{indexfig - (6 * (indexmontage - 1))} = imread([dirOutput(indexfig).folder '\' dirOutput(indexfig).name]);
                                    else
                                        I{indexfig - (6 * (indexmontage - 1))} = imread([dirOutput(indexfig).folder '/' dirOutput(indexfig).name]);
                                    end
                                end
                                
                                
                                
                                
                                % Creating fake images to fill the 3x2 image
                                for indexfig = (size(dirOutput,1) + 1) : (6 * indexmontage)
                                    I{indexfig - (6 * (indexmontage - 1))} = 255 * ones(size(I{1, 1}), 'uint8');
                                end
                                
                                Im = [I{1, 1} I{1, 2};
                                    I{1, 3} I{1, 4};
                                    I{1, 5} I{1, 6}];
                                imshow(Im, [], 'InitialMagnification','fit');
                                
                                if ispc
                                    % saveas(f, [handles.FIELDRT.env.userDir '\FIELDRT_Reports\Images\' module '\' struct_name_case{index1} '\' struct_name_temp '_montage_' int2str(indexmontage)] , 'tiff')
                                    % print(f, [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage)] , '-dtiffn', '-r500')
                                    imwrite(Im,[handles.FIELDRT.env.userDir '\FIELDRT_Reports\Images\' module '\' struct_name_case{index1} '\' struct_name_temp '_montage_' int2str(indexmontage + 1) '.tiff'],'tiff')
                                else
                                    % saveas(f, [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage)] , 'tiff')
                                    % print(f, [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage)] , '-dtiffn', '-r500')
                                    imwrite(Im,[handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage + 1) '.tiff'],'tiff')
                                end
                                
                                %             montage(fileNames, 'ThumbnailSize', [(size(I, 1) * 2) (size(I, 2) * 3)], 'Indices', (6 * (indexmontage - 1) + 1) : size(dirOutput,1), 'BackgroundColor', 'white', 'Size', [2 3]);
                                %
                                %             % Saving the montage
                                %             if ispc
                                %                 % saveas(f, [handles.FIELDRT.env.userDir '\FIELDRT_Reports\Images\' module '\' struct_name_case{index1} '\' struct_name_temp '_montage_' int2str(indexmontage)] , 'tiff')
                                %                 print(f, [handles.FIELDRT.env.userDir '\FIELDRT_Reports\Images\' module '\' struct_name_case{index1} '\' struct_name_temp '_montage_' int2str(indexmontage)] , '-dtiffn', '-r500')
                                %             else
                                %                 % saveas(f, [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage)] , 'tiff')
                                %                 print(f, [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage)] , '-dtiffn', '-r500')
                                %             end
                            else
                                % Only for the first structure a collage of four images is
                                % created. This is needed because of the layout of the first page of
                                % each functionality (i.e over/under, OARs...)
                                
                                
                                clear I
                                
                                for indexfig = (6 * (indexmontage - 1) + 1) : (6 * indexmontage)
                                    if ispc
                                        I{indexfig - (6 * (indexmontage - 1))} = imread([dirOutput(indexfig).folder '\' dirOutput(indexfig).name]);
                                    else
                                        I{indexfig - (6 * (indexmontage - 1))} = imread([dirOutput(indexfig).folder '/' dirOutput(indexfig).name]);
                                    end
                                end
                                
                                Im = [I{1, 1} I{1, 2};
                                    I{1, 3} I{1, 4};
                                    I{1, 5} I{1, 6}];
                                imshow(Im, [], 'InitialMagnification','fit');
                                
                                if ispc
                                    % saveas(f, [handles.FIELDRT.env.userDir '\FIELDRT_Reports\Images\' module '\' struct_name_case{index1} '\' struct_name_temp '_montage_' int2str(indexmontage)] , 'tiff')
                                    % print(f, [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage)] , '-dtiffn', '-r500')
                                    imwrite(Im,[handles.FIELDRT.env.userDir '\FIELDRT_Reports\Images\' module '\' struct_name_case{index1} '\' struct_name_temp '_montage_' int2str(indexmontage + 1) '.tiff'],'tiff')
                                else
                                    % saveas(f, [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage)] , 'tiff')
                                    % print(f, [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage)] , '-dtiffn', '-r500')
                                    imwrite(Im,[handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage + 1) '.tiff'],'tiff')
                                end
                                
                                
                                %             montage(fileNames, 'ThumbnailSize', [(size(I, 1) * 2) (size(I, 2) * 3)], 'Indices', (6 * (indexmontage - 1) + 1) : (6 * indexmontage), 'BackgroundColor', 'white', 'Size', [2 3]);
                                %
                                %             if ispc
                                %                 % saveas(f, [handles.FIELDRT.env.userDir '\FIELDRT_Reports\Images\' module '\' struct_name_case{index1} '\' struct_name_temp '_montage_' int2str(indexmontage)] , 'tiff')
                                %                 print(f, [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage)] , '-dtiffn', '-r500')
                                %             else
                                %                 % saveas(f, [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage)] , 'tiff')
                                %                 print(f, [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage)] , '-dtiffn', '-r500')
                                %             end
                            end
                        end
                    end
                    
                    indeximagestotal (1, index1) = nummont + 1;
                    
                elseif firstpageflag == 0
                    % creation of the collage for the other structures
                    if ispc
                        fileFolder = [handles.FIELDRT.env.userDir '\FIELDRT_Reports\Images\' module '\' struct_name_case{index1} '\'];
                    else
                        fileFolder = [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/'];
                    end
                    
                    dirOutput = dir(fullfile(fileFolder,[struct_name_case{index1} '*.tif']));
                    
                    % dirOutput = dirOutput(5 : end);
                    fileNames = string({dirOutput.name}); %retrieving filenames
                    
                    if size(dirOutput,1) > 0
                        % Numbers of montages --> we want 6 images per page
                        nummont = ceil(size(dirOutput,1) / 6);
                        
                        %     % read the first image to retrieve its resolution
                        %     I = imread(fileNames(1));
                        
                        % Creation of the collage with the other images
                        for indexmontage = 1 : nummont
                            
                            if indexmontage == nummont
                                
                                clear I
                                
                                % number of images not available to fill the 3*2 collage
                                % imadiff = (nummont * 6) - size(dirOutput,1);
                                
                                for indexfig = (6 * (indexmontage - 1) + 1) : size(dirOutput, 1)
                                    if ispc
                                        I{indexfig - (6 * (indexmontage - 1))} = imread([dirOutput(indexfig).folder '\' dirOutput(indexfig).name]);
                                    else
                                        I{indexfig - (6 * (indexmontage - 1))}  = imread([dirOutput(indexfig).folder '/' dirOutput(indexfig).name]);
                                    end
                                end
                                
                                
                                % Creating fake images to fill the 3x2 image
                                for indexfig = (size(dirOutput,1) + 1) : (6 * indexmontage)
                                    I{indexfig - (6 * (indexmontage - 1))} = 255 * ones(size(I{1, 1}), 'uint8');
                                end
                                
                                Im = [I{1, 1} I{1, 2};
                                    I{1, 3} I{1, 4};
                                    I{1, 5} I{1, 6}];
                                imshow(Im, [], 'InitialMagnification','fit');
                                
                                if ispc
                                    % saveas(f, [handles.FIELDRT.env.userDir '\FIELDRT_Reports\Images\' module '\' struct_name_case{index1} '\' struct_name_temp '_montage_' int2str(indexmontage)] , 'tiff')
                                    % print(f, [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage)] , '-dtiffn', '-r500')
                                    imwrite(Im,[handles.FIELDRT.env.userDir '\FIELDRT_Reports\Images\' module '\' struct_name_case{index1} '\' struct_name_temp '_montage_' int2str(indexmontage) '.tiff'],'tiff')
                                else
                                    % saveas(f, [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage)] , 'tiff')
                                    % print(f, [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage)] , '-dtiffn', '-r500')
                                    imwrite(Im,[handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage) '.tiff'],'tiff')
                                end
                                
                                %             montage(fileNames, 'ThumbnailSize', [(size(I, 1) * 2) (size(I, 2) * 3)], 'Indices', (6 * (indexmontage - 1) + 1) : size(dirOutput,1), 'BackgroundColor', 'white', 'Size', [2 3]);
                                %
                                %             % Saving the montage
                                %             if ispc
                                %                 % saveas(f, [handles.FIELDRT.env.userDir '\FIELDRT_Reports\Images\' module '\' struct_name_case{index1} '\' struct_name_temp '_montage_' int2str(indexmontage)] , 'tiff')
                                %                 print(f, [handles.FIELDRT.env.userDir '\FIELDRT_Reports\Images\' module '\' struct_name_case{index1} '\' struct_name_temp '_montage_' int2str(indexmontage)] , '-dtiffn', '-r500')
                                %             else
                                %                 % saveas(f, [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage)] , 'tiff')
                                %                 print(f, [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage)] , '-dtiffn', '-r500')
                                %             end
                            else
                                % Only for the first structure a collage of four images is
                                % created. This is needed because of the layout of the first page of
                                % each functionality (i.e over/under, OARs...)
                                
                                
                                clear I
                                
                                for indexfig = (6 * (indexmontage - 1) + 1) : (6 * indexmontage)
                                    if ispc
                                        I{indexfig - (6 * (indexmontage - 1))} = imread([dirOutput(indexfig).folder '\' dirOutput(indexfig).name]);
                                    else
                                        I{indexfig - (6 * (indexmontage - 1))}  = imread([dirOutput(indexfig).folder '/' dirOutput(indexfig).name]);
                                    end
                                end
                                
                                Im = [I{1, 1} I{1, 2};
                                    I{1, 3} I{1, 4};
                                    I{1, 5} I{1, 6}];
                                imshow(Im, [], 'InitialMagnification','fit');
                                
                                if ispc
                                    % saveas(f, [handles.FIELDRT.env.userDir '\FIELDRT_Reports\Images\' module '\' struct_name_case{index1} '\' struct_name_temp '_montage_' int2str(indexmontage)] , 'tiff')
                                    % print(f, [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage)] , '-dtiffn', '-r500')
                                    imwrite(Im,[handles.FIELDRT.env.userDir '\FIELDRT_Reports\Images\' module '\' struct_name_case{index1} '\' struct_name_temp '_montage_' int2str(indexmontage) '.tiff'],'tiff')
                                else
                                    % saveas(f, [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage)] , 'tiff')
                                    % print(f, [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage)] , '-dtiffn', '-r500')
                                    imwrite(Im,[handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage) '.tiff'],'tiff')
                                end
                                
                                
                                %             montage(fileNames, 'ThumbnailSize', [(size(I, 1) * 2) (size(I, 2) * 3)], 'Indices', (6 * (indexmontage - 1) + 1) : (6 * indexmontage), 'BackgroundColor', 'white', 'Size', [2 3]);
                                %
                                %             if ispc
                                %                 % saveas(f, [handles.FIELDRT.env.userDir '\FIELDRT_Reports\Images\' module '\' struct_name_case{index1} '\' struct_name_temp '_montage_' int2str(indexmontage)] , 'tiff')
                                %                 print(f, [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage)] , '-dtiffn', '-r500')
                                %             else
                                %                 % saveas(f, [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage)] , 'tiff')
                                %                 print(f, [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage)] , '-dtiffn', '-r500')
                                %             end
                            end
                        end
                    end
                    
                    % Number of figures for each structure
                    indeximagestotal (1, index1) = nummont;
                    
                    % Closing the figure
                    clear f
                end
            end
        case 'OARs'
            % Index for counting the images created
            indeximages_temp = 0;
            
            for index2 = range_total(index1, 1) : range_total(index1, 2)
                if module_imagetoshow(index1, index2) == 1 % The image was flagged so it needs to be shown
                    
                    hhl = [];
                    labelhhl = [];
                    
                    % Range of slices
                    minSlice = min(min(handles.planC{1, handles.scanfieldnum}.scanArray(:, :, index2)));
                    maxSlice = max(max(handles.planC{1, handles.scanfieldnum}.scanArray(:, :, index2)));
                    
                    % Retieving information for contours plotting
                    scale =  handles.planC{handles.scanfieldnum}.scanInfo(index2).grid1Units;
                    imageWidth(1, 1) =  handles.planC{handles.scanfieldnum}.scanInfo(index2).sizeOfDimension1;
                    imageWidth(2, 1) =  handles.planC{handles.scanfieldnum}.scanInfo(index2).sizeOfDimension2;
                    
                    if ~isempty(handles.planC{handles.scanfieldnum}.scanInfo(1).xOffset)
                        xCTOffset = handles.planC{handles.scanfieldnum}.scanInfo(1).xOffset;
                        yCTOffset = handles.planC{handles.scanfieldnum}.scanInfo(1).yOffset;
                    else
                        xCTOffset = 0;
                        yCTOffset = 0;
                    end
                    
                    pointsM = [];
                    
                    % Visualizing the image
                    hhh = imshow (handles.planC{1, handles.scanfieldnum}.scanArray(:, :, index2), [minSlice maxSlice], 'InitialMagnification','fit');
                    
                    hold on
                    
                    %% GS's structure
                    GSsstructure  = getStructNum([struct_name_temp '_GS'], handles.planC, indexS);
                    
                    hL1 = [];
                    
                    yCoordsGS = [];
                    xCoordsSG = [];
                    
                    if isempty(handles.planC{handles.structfieldnum}(GSsstructure).contour(index2).segments) 
                        
                        clear hL1
                        
                    else
                        
                        pointsM = handles.planC{handles.structfieldnum}(GSsstructure).contour(index2).segments.points;
                        
                        if isempty(pointsM) % Another check if there are points in the contour
                            
                        else
                            
                            xCoords = pointsM(:,1);
                            yCoords = pointsM(:,2);
                            
                            [rowV, colV] = aapmtom(xCoords/scale, yCoords/scale, xCTOffset/scale, yCTOffset/scale, imageWidth);
                            
                            rowV = [rowV; rowV(1)];
                            colV = [colV; colV(1)];
                            
                            % Plotting the contour
                            hL1 = plot(colV, rowV, 'Color', [1 1 0], 'linewidth', 4);
                        end
                    end
                    
                    % Legend for the contour
                    if exist('hL1')
                        if isempty(hL1)
                            
                        else
                            hhl = [hhl hL1];
                            labelhhl = [labelhhl {'Reference volume'}];
                        end
                    end
                    
                    %%  User's structure
                    USsstructure  = getStructNum(struct_name_temp, handles.planC, indexS);
                    
                    hL2 = [];
                    
                    yCoordsGS = [];
                    xCoordsSG = [];
                    
                    % Checking if there are points in the contour
                    if isempty(handles.planC{handles.structfieldnum}(USsstructure).contour(index2).segments) % Where there are no contours to plot
                        
                        clear hL2
                        
                    else
                        pointsM = handles.planC{handles.structfieldnum}(USsstructure).contour(index2).segments.points;
                        
                        if isempty(pointsM) % Another check if there are points in the contour
                            
                        else
                            xCoords = pointsM(:,1);
                            yCoords = pointsM(:,2);
                            
                            [rowV, colV] = aapmtom(xCoords/scale, yCoords/scale, xCTOffset/scale, yCTOffset/scale, imageWidth);
                            
                            rowV = [rowV; rowV(1)];
                            colV = [colV; colV(1)];
                            
                            % Plotting the contour
                            hL2 = plot(colV, rowV, 'Color', [0 1 1], 'linewidth', 4);
                        end
                    end
                    
                    % Legend for the contour
                    if exist('hL2')
                        if isempty(hL2)
                            
                        else
                            hhl = [hhl hL2];
                            labelhhl = [labelhhl, {'User volume'}];
                        end
                    else
                    end
                    
                    %% OARs (only flagged OARs will be shown)
                    % List of colours to use with OARs
                    colorlab = handles.FIELDRTGSCases(handles.Case(1,1)).OARs(handles.Case(1,2)).colorlab;
                    
                    flagOARstoshow = [];
                    indexflagOARstoshow = 1;
                    % To determine OARs to be shown
                    for indexoarmask = 1 : size(handles.OARsdummy_finalmask, 2)
                        tempmask_dummmyOARs = handles.OARsdummy_finalmask{index1, indexoarmask};
                        flag_label = 0; % flag to avoid adding in the legend two centroids derived from the same OAR
                        % Checking if there is an overlapping
                        CC_dummy(index1, indexoarmask) = bwconncomp(tempmask_dummmyOARs(:, :, index2)); % connectivity = 8 (default)
                        if CC_dummy(index1, indexoarmask).NumObjects > 0
                            % there is an intersection so the centroid is
                            % computed from the intersection between the US
                            % outline and the OAR (not the dummy OAR)
                            
                            % tempmask_OARs = handles.OARs_finalmask{index1, indexoarmask};
                            % CC(index1, indexoarmask) = bwconncomp(tempmask_OARs(:, :, index2)); % connectivity = 8 (default)
                            %
                            TT = regionprops(CC_dummy(index1, indexoarmask), 'Centroid', 'Area');
                            [val_TT, idx_TT] = max([TT.Area]);
                            
                            flagOARstoshow (indexflagOARstoshow, 1) = indexoarmask;   
                            indexflagOARstoshow = indexflagOARstoshow + 1;
                        else
                            % Nothing happens
                        end
                    end
                    
                    
                    for indexOARs = 1 : size(flagOARstoshow, 1) % Retrieving information about the selected structure
                        struct_name_temp_OARs = [struct_name_case_OARs{flagOARstoshow(indexOARs)} '_GS'];
                        OARsstructure  = getStructNum(struct_name_temp_OARs, handles.planC, indexS);
                        
                        hL4 = [];
                        
                        yCoordsGS = [];
                        xCoordsSG = [];
                        
                        % Checking if there are points in the contour
                        if (isempty(handles.planC{handles.structfieldnum}(OARsstructure).contour(index2).segments))  % Where there are no contours to plot
                            
                            clear hL4
                            
                        else
                            
                            pointsM = handles.planC{handles.structfieldnum}(OARsstructure).contour(index2).segments.points;
                            
                            if isempty(pointsM) % Another check if there are points in the contour
                                
                            else
                                
                                xCoords = pointsM(:,1);
                                yCoords = pointsM(:,2);
                                
                                [rowV, colV] = aapmtom(xCoords/scale, yCoords/scale, xCTOffset/scale, yCTOffset/scale, imageWidth);
                                
                                rowV = [rowV; rowV(1)];
                                colV = [colV; colV(1)];
                                
                                % Plotting the contour
                                hL4 = plot(colV, rowV, 'Color', colorlab(flagOARstoshow(indexOARs), :), 'linewidth', 4);
                            end
                        end
                        
                        % Legend for the contour
                        if exist('hL4')
                            if isempty(hL4)
                                
                            else
                                hhl = [hhl hL4];
                                %                                 struct_name_temp_OARs = struct_name_case_OARs{indexOARs};
                                %                                 labelhhl = [labelhhl {struct_name_temp_OARs(1 : end - 3)}];
                                labelhhl = [labelhhl struct_name_case_OARs{flagOARstoshow(indexOARs)}];
                            end
                        else
                        end
                    end
                    
                    
                    % Checking if there is an overlapping between the user selected
                    % structure and each OAR in the current slice. If yes, a flag appears.
                    clear CC
                    S = [];
                    
                    indexlegendline = 0;  % flag to avoid adding two blank lines in the legend more than one time
                    for indexoarmask = 1 : size(flagOARstoshow, 1)
                        tempmask_dummmyOARs = handles.OARsdummy_finalmask{index1, flagOARstoshow(indexoarmask)};
                        flag_label = 0; % flag to avoid adding in the legend two centroids derived from the same OAR
                        % Checking if there is an overlapping
                        CC_dummy(index1, indexoarmask) = bwconncomp(tempmask_dummmyOARs(:, :, index2)); % connectivity = 8 (default)
                        if CC_dummy(index1, indexoarmask).NumObjects > 0
                            % there is an intersection so the centroid is
                            % computed from the intersection between the US
                            % outline and the OAR (not the dummy OAR)
                            
                            % tempmask_OARs = handles.OARs_finalmask{index1, indexoarmask};
                            % CC(index1, indexoarmask) = bwconncomp(tempmask_OARs(:, :, index2)); % connectivity = 8 (default)
                            %
                            TT = regionprops(CC_dummy(index1, indexoarmask), 'Centroid', 'Area');
                            [val_TT, idx_TT] = max([TT.Area]);
                            
                            % To add two blank lines in the legend
                            if indexlegendline == 0
                                for indexx1 = 1 : 1 % Just one space added
                                    hL5 = plot([0.1 0.1], 'Color', 'none', 'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 1);
                                    hhl = [hhl, hL5];
                                    struct_name_temp_legend =  ' ';
                                    labelhhl = [labelhhl struct_name_temp_legend];
                                end
                                indexlegendline = 1;
                            end
                            
                            hLmask_new = plot(TT(idx_TT, 1).Centroid(1, 1), TT(idx_TT, 1).Centroid(1, 2), 'LineStyle', 'none', 'Marker', '.', 'color', colorlab(flagOARstoshow(indexoarmask), :),'MarkerSize', 40);
                            
                            if flag_label == 0 % so to avoid adding in the legend two centroids derived from the same OAR
                                hhl = [hhl, hLmask_new];
                                %                                 struct_name_temp_OARs = [struct_name_case_OARs{indexoarmask}];
                                %                                 labelhhl = [labelhhl , {['Too much ' struct_name_temp_OARs(1 : end - 3) ' included']}];
                                labelhhl = [labelhhl , {['Too much ' struct_name_case_OARs{flagOARstoshow(indexoarmask)} ' included']}];
                                flag_label = 1;
                            end
                        else
                            % Nothing appears
                        end
                    end
                    
                    if size(hhl, 2) >= 1
                        hhllegend = legend(hhl, labelhhl, 'FontSize', 32, 'TextColor', 'white', 'Color', 'Black');
                    else
                        
                    end
                    
                    % Title of the image
                    caption = sprintf('Slice %i', index2);
                    title(caption, 'FontSize', 40, 'Color' , 'k');
                    
                    indeximages_temp = indeximages_temp + 1;
                    
                    if ispc
                        % saveas(f, [handles.FIELDRT.env.userDir '\FIELDRT_Reports\Images\' module '\'  struct_name_case{index1} '\' struct_name_temp int2str(index2)] , 'tiff')
                        print(f, [handles.FIELDRT.env.userDir '\FIELDRT_Reports\Images\' module '\'  struct_name_case{index1} '\' struct_name_temp int2str(index2)] , '-dtiffn')
                    else
                        % saveas(f, [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp int2str(index2)] , 'tiff')
                        print(f, [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp int2str(index2)] , '-dtiffn')
                    end
                    
                    hold off
                end
            end
            
            % Index for counting the images created in each folder
            % indeximages_total(2) = indeximages_temp;
            
            % closing the figure
            clear f
            
            
            %% Creating the montage with the created figures
            
            %     % figure for the montage
            %     f = figure('visible','off', 'Position', [screen(1) screen(2) (screen(3) + 500) (screen(4) + 200) ]);
            %
            %     % To have the figure with Minimal White Space
            %     ax = gca;
            %     outerpos = ax.OuterPosition;
            %     ti = ax.TightInset;
            %     left = outerpos(1) + ti(1);
            %     bottom = outerpos(2) + ti(2);
            %     ax_width = outerpos(3) - ti(1) - ti(3);
            %     ax_height = outerpos(4) - ti(2) - ti(4);
            %     ax.Position = [left bottom ax_width ax_height];

            
            if indeximages_temp == 0 % no images are created so no montage to create and save
                indeximagestotal (1, index1) = 0;
            else % there are images so a montage need to be created
                
                nummont = 0;
                
                if firstpageflag == 1 %Only for the first page of each section there are 4 images. The other pages all have 6 images.
                    if ispc
                        fileFolder = [handles.FIELDRT.env.userDir '\FIELDRT_Reports\Images\' module '\' struct_name_case{index1} '\'];
                    else
                        fileFolder = [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/'];
                    end
                    
                    dirOutput_iniz = dir(fullfile(fileFolder,[struct_name_case{index1} '*.tif']));
                    
                    % Only for the first structure a collage of four images is
                    % created. This is needed because of the layout of the first page of
                    % each functionality (i.e over/under, OARs...)
                    if size(dirOutput_iniz, 1) < 4
                        dirOutput_temp = dirOutput_iniz(1 : size(dirOutput_iniz, 1));
                        fileNames_temp = string({dirOutput_temp.name}); %retrieving filenames
                        firstpageflag = 0;
                    else
                        dirOutput_temp = dirOutput_iniz(1 : 4);
                        fileNames_temp = string({dirOutput_temp.name}); %retrieving filenames
                        dirOutput = dirOutput_iniz(5 : end);
                        firstpageflag = 0;
                    end
                    
                    % Creation of the first collage with four images for the first
                    % structure
                    clear I
                    
                    if size(dirOutput_temp, 1) > 0
                        for indexfig = 1 : size(dirOutput_temp, 1)
                            if ispc
                                I{indexfig} = imread([dirOutput_temp(indexfig).folder '\' dirOutput_temp(indexfig).name]);
                            else
                                I{indexfig} = imread([dirOutput_temp(indexfig).folder '/' dirOutput_temp(indexfig).name]);
                            end
                        end
                        
                        if size(dirOutput_temp, 1) == 1
                             Im = [I{1, 1} ones(size(I{1, 1}))*255;
                                 ones(size(I{1, 1}))*255 ones(size(I{1, 1}))*255];
                        end
                        
                        if size(dirOutput_temp, 1) == 2
                             Im = [I{1, 1} I{1, 2};
                                 ones(size(I{1, 1}))*255 ones(size(I{1, 1}))*255];
                        end
                        
                        if size(dirOutput_temp, 1) == 3
                             Im = [I{1, 1} I{1, 2};
                            I{1, 3} ones(size(I{1, 1}))*255];
                        end
                        
                        if size(dirOutput_temp, 1) == 4
                             Im = [I{1, 1} I{1, 2};
                            I{1, 3} I{1, 4}];
                        end
                       
                        % imshow(Im, [], 'InitialMagnification','fit');
                        
                        if ispc
                            % saveas(f, [handles.FIELDRT.env.userDir '\FIELDRT_Reports\Images\' module '\' struct_name_case{index1} '\' struct_name_temp '_montage_' int2str(indexmontage)] , 'tiff')
                            % print(f, [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage)] , '-dtiffn', '-r500')
                            imwrite(Im,[handles.FIELDRT.env.userDir '\FIELDRT_Reports\Images\' module '\' struct_name_case{index1} '\' struct_name_temp '_montage_' int2str(1) '.tiff'],'tiff')
                        else
                            % saveas(f, [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage)] , 'tiff')
                            % print(f, [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage)] , '-dtiffn', '-r500')
                            imwrite(Im,[handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(1) '.tiff'],'tiff')
                        end
                    end
                    
                    
                    % creation of the montage with 6 images
                    if size(dirOutput_iniz,1) > 4 
                        
                        nummont = ceil(size(dirOutput,1) / 6);
                        
                        for indexmontage = 1 : nummont
                            
                            if indexmontage == nummont
                                
                                clear I
                                
                                % number of images not available to fill the 3*2 collage
                                % imadiff = (nummont * 6) - size(dirOutput,1);
                                
                                for indexfig = (6 * (indexmontage - 1) + 1) : size(dirOutput, 1)
                                    if ispc
                                        I{indexfig - (6 * (indexmontage - 1))} = imread([dirOutput(indexfig).folder '\' dirOutput(indexfig).name]);
                                    else
                                        I{indexfig - (6 * (indexmontage - 1))} = imread([dirOutput(indexfig).folder '/' dirOutput(indexfig).name]);
                                    end
                                end
                                
                                
                                
                                
                                % Creating fake images to fill the 3x2 image
                                for indexfig = (size(dirOutput,1) + 1) : (6 * indexmontage)
                                    I{indexfig - (6 * (indexmontage - 1))} = 255 * ones(size(I{1, 1}), 'uint8');
                                end
                                
                                Im = [I{1, 1} I{1, 2};
                                    I{1, 3} I{1, 4};
                                    I{1, 5} I{1, 6}];
                                imshow(Im, [], 'InitialMagnification','fit');
                                
                                if ispc
                                    % saveas(f, [handles.FIELDRT.env.userDir '\FIELDRT_Reports\Images\' module '\' struct_name_case{index1} '\' struct_name_temp '_montage_' int2str(indexmontage)] , 'tiff')
                                    % print(f, [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage)] , '-dtiffn', '-r500')
                                    imwrite(Im,[handles.FIELDRT.env.userDir '\FIELDRT_Reports\Images\' module '\' struct_name_case{index1} '\' struct_name_temp '_montage_' int2str(indexmontage + 1) '.tiff'],'tiff')
                                else
                                    % saveas(f, [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage)] , 'tiff')
                                    % print(f, [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage)] , '-dtiffn', '-r500')
                                    imwrite(Im,[handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage + 1) '.tiff'],'tiff')
                                end
                                
                                %             montage(fileNames, 'ThumbnailSize', [(size(I, 1) * 2) (size(I, 2) * 3)], 'Indices', (6 * (indexmontage - 1) + 1) : size(dirOutput,1), 'BackgroundColor', 'white', 'Size', [2 3]);
                                %
                                %             % Saving the montage
                                %             if ispc
                                %                 % saveas(f, [handles.FIELDRT.env.userDir '\FIELDRT_Reports\Images\' module '\' struct_name_case{index1} '\' struct_name_temp '_montage_' int2str(indexmontage)] , 'tiff')
                                %                 print(f, [handles.FIELDRT.env.userDir '\FIELDRT_Reports\Images\' module '\' struct_name_case{index1} '\' struct_name_temp '_montage_' int2str(indexmontage)] , '-dtiffn', '-r500')
                                %             else
                                %                 % saveas(f, [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage)] , 'tiff')
                                %                 print(f, [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage)] , '-dtiffn', '-r500')
                                %             end
                            else
                                % Only for the first structure a collage of four images is
                                % created. This is needed because of the layout of the first page of
                                % each functionality (i.e over/under, OARs...)
                                
                                
                                clear I
                                
                                for indexfig = (6 * (indexmontage - 1) + 1) : (6 * indexmontage)
                                    if ispc
                                        I{indexfig - (6 * (indexmontage - 1))} = imread([dirOutput(indexfig).folder '\' dirOutput(indexfig).name]);
                                    else
                                        I{indexfig - (6 * (indexmontage - 1))} = imread([dirOutput(indexfig).folder '/' dirOutput(indexfig).name]);
                                    end
                                end
                                
                                Im = [I{1, 1} I{1, 2};
                                    I{1, 3} I{1, 4};
                                    I{1, 5} I{1, 6}];
                                imshow(Im, [], 'InitialMagnification','fit');
                                
                                if ispc
                                    % saveas(f, [handles.FIELDRT.env.userDir '\FIELDRT_Reports\Images\' module '\' struct_name_case{index1} '\' struct_name_temp '_montage_' int2str(indexmontage)] , 'tiff')
                                    % print(f, [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage)] , '-dtiffn', '-r500')
                                    imwrite(Im,[handles.FIELDRT.env.userDir '\FIELDRT_Reports\Images\' module '\' struct_name_case{index1} '\' struct_name_temp '_montage_' int2str(indexmontage + 1) '.tiff'],'tiff')
                                else
                                    % saveas(f, [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage)] , 'tiff')
                                    % print(f, [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage)] , '-dtiffn', '-r500')
                                    imwrite(Im,[handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage + 1) '.tiff'],'tiff')
                                end
                                
                                
                                %             montage(fileNames, 'ThumbnailSize', [(size(I, 1) * 2) (size(I, 2) * 3)], 'Indices', (6 * (indexmontage - 1) + 1) : (6 * indexmontage), 'BackgroundColor', 'white', 'Size', [2 3]);
                                %
                                %             if ispc
                                %                 % saveas(f, [handles.FIELDRT.env.userDir '\FIELDRT_Reports\Images\' module '\' struct_name_case{index1} '\' struct_name_temp '_montage_' int2str(indexmontage)] , 'tiff')
                                %                 print(f, [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage)] , '-dtiffn', '-r500')
                                %             else
                                %                 % saveas(f, [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage)] , 'tiff')
                                %                 print(f, [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage)] , '-dtiffn', '-r500')
                                %             end
                            end
                        end
                    end
                    
                    indeximagestotal (1, index1) = nummont + 1;
                    
                elseif firstpageflag == 0 
                    % creation of the collage for the other structures
                    if ispc
                        fileFolder = [handles.FIELDRT.env.userDir '\FIELDRT_Reports\Images\' module '\' struct_name_case{index1} '\'];
                    else
                        fileFolder = [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/'];
                    end
                    
                    dirOutput_iniz = dir(fullfile(fileFolder,[struct_name_case{index1} '*.tif']));
                    
                                      
                    
                    dirOutput = dir(fullfile(fileFolder,[struct_name_case{index1} '*.tif']));
                    
                    % dirOutput = dirOutput(5 : end);
                    fileNames = string({dirOutput.name}); %retrieving filenames
                    
                    if size(dirOutput,1) > 0
                        % Numbers of montages --> we want 6 images per page
                        nummont = ceil(size(dirOutput,1) / 6);
                        
                        %     % read the first image to retrieve its resolution
                        %     I = imread(fileNames(1));
                        
                        % Creation of the collage with the other images
                        for indexmontage = 1 : nummont
                            
                            if indexmontage == nummont
                                
                                clear I
                                
                                % number of images not available to fill the 3*2 collage
                                % imadiff = (nummont * 6) - size(dirOutput,1);
                                
                                for indexfig = (6 * (indexmontage - 1) + 1) : size(dirOutput, 1)
                                    if ispc
                                        I{indexfig - (6 * (indexmontage - 1))} = imread([dirOutput(indexfig).folder '\' dirOutput(indexfig).name]);
                                    else
                                        I{indexfig - (6 * (indexmontage - 1))}  = imread([dirOutput(indexfig).folder '/' dirOutput(indexfig).name]);
                                    end
                                end
                                
                                
                                % Creating fake images to fill the 3x2 image
                                for indexfig = (size(dirOutput,1) + 1) : (6 * indexmontage)
                                    I{indexfig - (6 * (indexmontage - 1))} = 255 * ones(size(I{1, 1}), 'uint8');
                                end
                                
                                Im = [I{1, 1} I{1, 2};
                                    I{1, 3} I{1, 4};
                                    I{1, 5} I{1, 6}];
                                imshow(Im, [], 'InitialMagnification','fit');
                                
                                if ispc
                                    % saveas(f, [handles.FIELDRT.env.userDir '\FIELDRT_Reports\Images\' module '\' struct_name_case{index1} '\' struct_name_temp '_montage_' int2str(indexmontage)] , 'tiff')
                                    % print(f, [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage)] , '-dtiffn', '-r500')
                                    imwrite(Im,[handles.FIELDRT.env.userDir '\FIELDRT_Reports\Images\' module '\' struct_name_case{index1} '\' struct_name_temp '_montage_' int2str(indexmontage) '.tiff'],'tiff')
                                else
                                    % saveas(f, [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage)] , 'tiff')
                                    % print(f, [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage)] , '-dtiffn', '-r500')
                                    imwrite(Im,[handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage) '.tiff'],'tiff')
                                end
                                
                                %             montage(fileNames, 'ThumbnailSize', [(size(I, 1) * 2) (size(I, 2) * 3)], 'Indices', (6 * (indexmontage - 1) + 1) : size(dirOutput,1), 'BackgroundColor', 'white', 'Size', [2 3]);
                                %
                                %             % Saving the montage
                                %             if ispc
                                %                 % saveas(f, [handles.FIELDRT.env.userDir '\FIELDRT_Reports\Images\' module '\' struct_name_case{index1} '\' struct_name_temp '_montage_' int2str(indexmontage)] , 'tiff')
                                %                 print(f, [handles.FIELDRT.env.userDir '\FIELDRT_Reports\Images\' module '\' struct_name_case{index1} '\' struct_name_temp '_montage_' int2str(indexmontage)] , '-dtiffn', '-r500')
                                %             else
                                %                 % saveas(f, [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage)] , 'tiff')
                                %                 print(f, [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage)] , '-dtiffn', '-r500')
                                %             end
                            else
                                % Only for the first structure a collage of four images is
                                % created. This is needed because of the layout of the first page of
                                % each functionality (i.e over/under, OARs...)
                                
                                
                                clear I
                                
                                for indexfig = (6 * (indexmontage - 1) + 1) : (6 * indexmontage)
                                    if ispc
                                        I{indexfig - (6 * (indexmontage - 1))} = imread([dirOutput(indexfig).folder '\' dirOutput(indexfig).name]);
                                    else
                                        I{indexfig - (6 * (indexmontage - 1))}  = imread([dirOutput(indexfig).folder '/' dirOutput(indexfig).name]);
                                    end
                                end
                                
                                Im = [I{1, 1} I{1, 2};
                                    I{1, 3} I{1, 4};
                                    I{1, 5} I{1, 6}];
                                imshow(Im, [], 'InitialMagnification','fit');
                                
                                if ispc
                                    % saveas(f, [handles.FIELDRT.env.userDir '\FIELDRT_Reports\Images\' module '\' struct_name_case{index1} '\' struct_name_temp '_montage_' int2str(indexmontage)] , 'tiff')
                                    % print(f, [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage)] , '-dtiffn', '-r500')
                                    imwrite(Im,[handles.FIELDRT.env.userDir '\FIELDRT_Reports\Images\' module '\' struct_name_case{index1} '\' struct_name_temp '_montage_' int2str(indexmontage) '.tiff'],'tiff')
                                else
                                    % saveas(f, [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage)] , 'tiff')
                                    % print(f, [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage)] , '-dtiffn', '-r500')
                                    imwrite(Im,[handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage) '.tiff'],'tiff')
                                end
                                
                                
                                %             montage(fileNames, 'ThumbnailSize', [(size(I, 1) * 2) (size(I, 2) * 3)], 'Indices', (6 * (indexmontage - 1) + 1) : (6 * indexmontage), 'BackgroundColor', 'white', 'Size', [2 3]);
                                %
                                %             if ispc
                                %                 % saveas(f, [handles.FIELDRT.env.userDir '\FIELDRT_Reports\Images\' module '\' struct_name_case{index1} '\' struct_name_temp '_montage_' int2str(indexmontage)] , 'tiff')
                                %                 print(f, [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage)] , '-dtiffn', '-r500')
                                %             else
                                %                 % saveas(f, [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage)] , 'tiff')
                                %                 print(f, [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage)] , '-dtiffn', '-r500')
                                %             end
                            end
                        end
                    end
                    
                    % Number of figures for each structure
                    indeximagestotal (1, index1) = nummont;
                    
                    % Closing the figure
                    clear f
                end
            end
        case 'Minmax'
            % Index for counting the images created
            indeximages_temp = 0;
            
            for index2 = range_total(index1, 1) : range_total(index1, 2)
                
                if module_imagetoshow(index1, index2) == 1 % The image was flagged so it needs to be shown
                    hhl = [];
                    labelhhl = [];
                    
                    % Range of slices
                    minSlice = min(min(handles.planC{1, handles.scanfieldnum}.scanArray(:, :, index2)));
                    maxSlice = max(max(handles.planC{1, handles.scanfieldnum}.scanArray(:, :, index2)));
                    
                    % Retieving information for contours plotting
                    scale =  handles.planC{handles.scanfieldnum}.scanInfo(index2).grid1Units;
                    imageWidth(1, 1) =  handles.planC{handles.scanfieldnum}.scanInfo(index2).sizeOfDimension1;
                    imageWidth(2, 1) =  handles.planC{handles.scanfieldnum}.scanInfo(index2).sizeOfDimension2;
                    
                    if ~isempty(handles.planC{handles.scanfieldnum}.scanInfo(1).xOffset)
                        xCTOffset = handles.planC{handles.scanfieldnum}.scanInfo(1).xOffset;
                        yCTOffset = handles.planC{handles.scanfieldnum}.scanInfo(1).yOffset;
                    else
                        xCTOffset = 0;
                        yCTOffset = 0;
                    end
                    
                    pointsM = [];
                    
                    % Visualizing the image
                    hhh = imshow (handles.planC{1, handles.scanfieldnum}.scanArray(:, :, index2), [minSlice maxSlice], 'InitialMagnification','fit');
                    
                    hold on
                    
                    %% GS's structure
                    GSsstructure  = getStructNum([struct_name_temp '_GS'], handles.planC, indexS);
                    
                    hL1 = [];
                    
                    yCoordsGS = [];
                    xCoordsSG = [];
                    
                    % Checking if there are points in the contour
                    if (isempty(handles.planC{handles.structfieldnum}(GSsstructure).contour(index2).segments))  % Where there are no contours to plot
                        
                        clear hL1
                    else
                        
                        pointsM = handles.planC{handles.structfieldnum}(GSsstructure).contour(index2).segments.points;
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
                            
                            hL1 = plot(colV, rowV, 'Color', [1 1 0], 'linewidth', 4);
                        end
                    end
                    
                    % Legend for the contour
                    if exist('hL1')
                        if isempty(hL1)
                            
                        else
                            hhl = [hhl hL1];
                            labelhhl = [labelhhl {'Reference volume'}];
                        end
                    else
                        
                    end
                    
                    
                    % User's structure
                    
                    USsstructure  = getStructNum(struct_name_temp, handles.planC, indexS);
                    
                    hL2 = [];
                    
                    yCoordsGS = [];
                    xCoordsSG = [];
                    
                    if isempty(handles.planC{handles.structfieldnum}(USsstructure).contour(index2).segments)
                        xCoords = [];
                        yCoords = [];
                        
                        hhl = hhl;
                        
                        clear hL2 % Making sure there is no vector of chart line objects
                    else
                        pointsM = handles.planC{handles.structfieldnum}(USsstructure).contour(index2).segments.points;
                        
                        if isempty(pointsM) % Where there are no contours to plot
                            xCoords = [];
                            yCoords = [];
                            
                            hhl = hhl;
                            
                        else
                            xCoords = pointsM(:,1);
                            yCoords = pointsM(:,2);
                            
                            [rowV, colV] = aapmtom(xCoords/scale, yCoords/scale, xCTOffset/scale, yCTOffset/scale, imageWidth);
                            
                            rowV = [rowV; rowV(1)];
                            colV = [colV; colV(1)];
                            
                            % Plotting the contour
                            hL2 = plot(colV, rowV, 'Color', [0 1 1], 'linewidth', 4);
                        end
                    end
                    
                    % Legend for the contour
                    if exist('hL2')
                        if isempty(hL2)
                            
                        else
                            hhl = [hhl hL2];
                            labelhhl = [labelhhl {'User volume'}];
                        end
                    else
                        
                    end
                    
                    
                    % GS's structure max
                    GSsstructuremax  = getStructNum([struct_name_temp '_max_GS'], handles.planC, indexS);
                    
                    if (isempty(handles.planC{handles.structfieldnum}(GSsstructuremax).contour(index2).segments))  % Where there are no contours to plot
                        
                        clear hL4 % Making sure there is no vector of chart line objects
                        
                    else
                        
                        pointsM = handles.planC{handles.structfieldnum}(GSsstructuremax).contour(index2).segments.points;
                        
                        if isempty(pointsM) % Where there are no contours to plot
                            xCoords = [];
                            yCoords = [];
                            
                            hhl = hhl;
                            
                        else
                            xCoords = pointsM(:,1);
                            yCoords = pointsM(:,2);
                            
                            [rowV, colV] = aapmtom(xCoords/scale, yCoords/scale, xCTOffset/scale, yCTOffset/scale, imageWidth);
                            
                            rowV = [rowV; rowV(1)];
                            colV = [colV; colV(1)];
                            
                            % Plotting the contour
                            hL4 = plot(colV, rowV, 'Color', [1 0 0], 'linewidth', 4);
                        end
                    end
                    
                    % Legend for the contour
                    if exist('hL4')
                        if isempty(hL4)
                            
                        else
                            hhl = [hhl hL4];
                            labelhhl = [labelhhl {[struct_name_temp ' max']}];
                        end
                    else
                        
                    end
                    
                    % GS's structure min
                    GSsstructuremin  = getStructNum([struct_name_temp '_min_GS'], handles.planC, indexS);
                    
                    if (isempty(handles.planC{handles.structfieldnum}(GSsstructuremin).contour(index2).segments))  % Where there are no contours to plot
                        
                        clear hL5 % Making sure there is no vector of chart line objects
                    else
                        
                        pointsM = handles.planC{handles.structfieldnum}(GSsstructuremin).contour(index2).segments.points;
                        
                        if isempty(pointsM) % Where there are no contours to plot
                            xCoords = [];
                            yCoords = [];
                            
                            hhl = hhl;
                            
                        else
                            xCoords = pointsM(:,1);
                            yCoords = pointsM(:,2);
                            
                            [rowV, colV] = aapmtom(xCoords/scale, yCoords/scale, xCTOffset/scale, yCTOffset/scale, imageWidth);
                            
                            rowV = [rowV; rowV(1)];
                            colV = [colV; colV(1)];
                            
                            % Plotting the contour
                            hL5 = plot(colV, rowV, 'Color', [0 1 0], 'linewidth', 4);
                        end
                    end
                    
                    % Legend for the contour
                    if exist('hL5')
                        if isempty(hL5)
                            
                        else
                            hhl = [hhl hL5];
                            labelhhl = [labelhhl {[struct_name_temp ' min']}];
                        end
                    else
                        
                    end
                    
                    if size(hhl, 2) >= 1
                        hhllegend = legend(hhl, labelhhl, 'FontSize', 32, 'TextColor', 'white', 'Color', 'Black');
                    else
                        
                    end
                    
                    % Title of the image
                    caption = sprintf('Slice %i', index2);
                    title(caption, 'FontSize', 40, 'Color' , 'k');
                    
                    indeximages_temp = indeximages_temp + 1;
                    
                    if ispc
                        % saveas(f, [handles.FIELDRT.env.userDir '\FIELDRT_Reports\Images\' module '\'  struct_name_case{index1} '\' struct_name_temp int2str(index2)] , 'tiff')
                        print(f, [handles.FIELDRT.env.userDir '\FIELDRT_Reports\Images\' module '\'  struct_name_case{index1} '\' struct_name_temp int2str(index2)] , '-dtiffn')
                    else
                        % saveas(f, [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp int2str(index2)] , 'tiff')
                        print(f, [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp int2str(index2)] , '-dtiffn')
                    end
                    
                    hold off
                end
            end
            
            % Index for counting the images created in each folder
            % indeximages_total(3) = indeximages_temp;
            
            %% Creating the montage with the created figures
            
            %     % figure for the montage
            %     f = figure('visible','off', 'Position', [screen(1) screen(2) (screen(3) + 500) (screen(4) + 200) ]);
            %
            %     % To have the figure with Minimal White Space
            %     ax = gca;
            %     outerpos = ax.OuterPosition;
            %     ti = ax.TightInset;
            %     left = outerpos(1) + ti(1);
            %     bottom = outerpos(2) + ti(2);
            %     ax_width = outerpos(3) - ti(1) - ti(3);
            %     ax_height = outerpos(4) - ti(2) - ti(4);
            %     ax.Position = [left bottom ax_width ax_height];
            
            if indeximages_temp == 0 % no images are created so no montage to create and save
                indeximagestotal (1, index1) = 0;
            else % there are images so a montage need to be created
                
                nummont = 0;
                
                if firstpageflag == 1
                    if ispc
                        fileFolder = [handles.FIELDRT.env.userDir '\FIELDRT_Reports\Images\' module '\' struct_name_case{index1} '\'];
                    else
                        fileFolder = [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/'];
                    end
                    
                    dirOutput_iniz = dir(fullfile(fileFolder,[struct_name_case{index1} '*.tif']));
                    
                    % Only for the first structure a collage of four images is
                    % created. This is needed because of the layout of the first page of
                    % each functionality (i.e over/under, OARs...)
                    if size(dirOutput_iniz, 1) < 4
                        dirOutput_temp = dirOutput_iniz(1 : size(dirOutput_iniz, 1));
                        fileNames_temp = string({dirOutput_temp.name}); %retrieving filenames
                        firstpageflag = 0;
                    else
                        dirOutput_temp = dirOutput_iniz(1 : 4);
                        fileNames_temp = string({dirOutput_temp.name}); %retrieving filenames
                        dirOutput = dirOutput_iniz(5 : end);
                        firstpageflag = 0;
                    end                    
                    
                    % Creation of the first collage with four images for the first
                    % structure
                    clear I
                    
                    if size(dirOutput_temp, 1) > 0
                        for indexfig = 1 : size(dirOutput_temp, 1)
                            if ispc
                                I{indexfig} = imread([dirOutput_temp(indexfig).folder '\' dirOutput_temp(indexfig).name]);
                            else
                                I{indexfig} = imread([dirOutput_temp(indexfig).folder '/' dirOutput_temp(indexfig).name]);
                            end
                        end
                        
                        if size(dirOutput_temp, 1) == 1
                             Im = [I{1, 1} ones(size(I{1, 1}))*255;
                                 ones(size(I{1, 1}))*255 ones(size(I{1, 1}))*255];
                        end
                        
                        if size(dirOutput_temp, 1) == 2
                             Im = [I{1, 1} I{1, 2};
                                 ones(size(I{1, 1}))*255 ones(size(I{1, 1}))*255];
                        end
                        
                        if size(dirOutput_temp, 1) == 3
                             Im = [I{1, 1} I{1, 2};
                            I{1, 3} ones(size(I{1, 1}))*255];
                        end
                        
                        if size(dirOutput_temp, 1) == 4
                             Im = [I{1, 1} I{1, 2};
                            I{1, 3} I{1, 4}];
                        end
                        
                        if ispc
                            % saveas(f, [handles.FIELDRT.env.userDir '\FIELDRT_Reports\Images\' module '\' struct_name_case{index1} '\' struct_name_temp '_montage_' int2str(indexmontage)] , 'tiff')
                            % print(f, [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage)] , '-dtiffn', '-r500')
                            imwrite(Im,[handles.FIELDRT.env.userDir '\FIELDRT_Reports\Images\' module '\' struct_name_case{index1} '\' struct_name_temp '_montage_' int2str(1) '.tiff'],'tiff')
                        else
                            % saveas(f, [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage)] , 'tiff')
                            % print(f, [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage)] , '-dtiffn', '-r500')
                            imwrite(Im,[handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(1) '.tiff'],'tiff')
                        end
                    end
                    
                    
                    % creation of the montage with 6 images
                    if size(dirOutput,1) > 0
                        
                        nummont = ceil(size(dirOutput,1) / 6);
                        
                        for indexmontage = 1 : nummont
                            
                            if indexmontage == nummont
                                
                                clear I
                                
                                % number of images not available to fill the 3*2 collage
                                % imadiff = (nummont * 6) - size(dirOutput,1);
                                
                                for indexfig = (6 * (indexmontage - 1) + 1) : size(dirOutput, 1)
                                    if ispc
                                        I{indexfig - (6 * (indexmontage - 1))} = imread([dirOutput(indexfig).folder '\' dirOutput(indexfig).name]);
                                    else
                                        I{indexfig - (6 * (indexmontage - 1))} = imread([dirOutput(indexfig).folder '/' dirOutput(indexfig).name]);
                                    end
                                end
                                
                                
                                
                                
                                % Creating fake images to fill the 3x2 image
                                for indexfig = (size(dirOutput,1) + 1) : (6 * indexmontage)
                                    I{indexfig - (6 * (indexmontage - 1))} = 255 * ones(size(I{1, 1}), 'uint8');
                                end
                                
                                Im = [I{1, 1} I{1, 2};
                                    I{1, 3} I{1, 4};
                                    I{1, 5} I{1, 6}];
                                imshow(Im, [], 'InitialMagnification','fit');
                                
                                if ispc
                                    % saveas(f, [handles.FIELDRT.env.userDir '\FIELDRT_Reports\Images\' module '\' struct_name_case{index1} '\' struct_name_temp '_montage_' int2str(indexmontage)] , 'tiff')
                                    % print(f, [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage)] , '-dtiffn', '-r500')
                                    imwrite(Im,[handles.FIELDRT.env.userDir '\FIELDRT_Reports\Images\' module '\' struct_name_case{index1} '\' struct_name_temp '_montage_' int2str(indexmontage + 1) '.tiff'],'tiff')
                                else
                                    % saveas(f, [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage)] , 'tiff')
                                    % print(f, [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage)] , '-dtiffn', '-r500')
                                    imwrite(Im,[handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage + 1) '.tiff'],'tiff')
                                end
                                
                                %             montage(fileNames, 'ThumbnailSize', [(size(I, 1) * 2) (size(I, 2) * 3)], 'Indices', (6 * (indexmontage - 1) + 1) : size(dirOutput,1), 'BackgroundColor', 'white', 'Size', [2 3]);
                                %
                                %             % Saving the montage
                                %             if ispc
                                %                 % saveas(f, [handles.FIELDRT.env.userDir '\FIELDRT_Reports\Images\' module '\' struct_name_case{index1} '\' struct_name_temp '_montage_' int2str(indexmontage)] , 'tiff')
                                %                 print(f, [handles.FIELDRT.env.userDir '\FIELDRT_Reports\Images\' module '\' struct_name_case{index1} '\' struct_name_temp '_montage_' int2str(indexmontage)] , '-dtiffn', '-r500')
                                %             else
                                %                 % saveas(f, [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage)] , 'tiff')
                                %                 print(f, [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage)] , '-dtiffn', '-r500')
                                %             end
                            else
                                % Only for the first structure a collage of four images is
                                % created. This is needed because of the layout of the first page of
                                % each functionality (i.e over/under, OARs...)
                                
                                
                                clear I
                                
                                for indexfig = (6 * (indexmontage - 1) + 1) : (6 * indexmontage)
                                    if ispc
                                        I{indexfig - (6 * (indexmontage - 1))} = imread([dirOutput(indexfig).folder '\' dirOutput(indexfig).name]);
                                    else
                                        I{indexfig - (6 * (indexmontage - 1))} = imread([dirOutput(indexfig).folder '/' dirOutput(indexfig).name]);
                                    end
                                end
                                
                                Im = [I{1, 1} I{1, 2};
                                    I{1, 3} I{1, 4};
                                    I{1, 5} I{1, 6}];
                                imshow(Im, [], 'InitialMagnification','fit');
                                
                                if ispc
                                    % saveas(f, [handles.FIELDRT.env.userDir '\FIELDRT_Reports\Images\' module '\' struct_name_case{index1} '\' struct_name_temp '_montage_' int2str(indexmontage)] , 'tiff')
                                    % print(f, [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage)] , '-dtiffn', '-r500')
                                    imwrite(Im,[handles.FIELDRT.env.userDir '\FIELDRT_Reports\Images\' module '\' struct_name_case{index1} '\' struct_name_temp '_montage_' int2str(indexmontage + 1) '.tiff'],'tiff')
                                else
                                    % saveas(f, [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage)] , 'tiff')
                                    % print(f, [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage)] , '-dtiffn', '-r500')
                                    imwrite(Im,[handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage + 1) '.tiff'],'tiff')
                                end
                                
                                
                                %             montage(fileNames, 'ThumbnailSize', [(size(I, 1) * 2) (size(I, 2) * 3)], 'Indices', (6 * (indexmontage - 1) + 1) : (6 * indexmontage), 'BackgroundColor', 'white', 'Size', [2 3]);
                                %
                                %             if ispc
                                %                 % saveas(f, [handles.FIELDRT.env.userDir '\FIELDRT_Reports\Images\' module '\' struct_name_case{index1} '\' struct_name_temp '_montage_' int2str(indexmontage)] , 'tiff')
                                %                 print(f, [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage)] , '-dtiffn', '-r500')
                                %             else
                                %                 % saveas(f, [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage)] , 'tiff')
                                %                 print(f, [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage)] , '-dtiffn', '-r500')
                                %             end
                            end
                        end
                    end
                    
                    indeximagestotal (1, index1) = nummont + 1;
                    
                elseif firstpageflag == 0 
                    % creation of the collage for the other structures
                    if ispc
                        fileFolder = [handles.FIELDRT.env.userDir '\FIELDRT_Reports\Images\' module '\' struct_name_case{index1} '\'];
                    else
                        fileFolder = [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/'];
                    end
                    
                    dirOutput = dir(fullfile(fileFolder,[struct_name_case{index1} '*.tif']));
                    
                    % dirOutput = dirOutput(5 : end);
                    fileNames = string({dirOutput.name}); %retrieving filenames
                    
                    if size(dirOutput,1) > 0
                        % Numbers of montages --> we want 6 images per page
                        nummont = ceil(size(dirOutput,1) / 6);
                        
                        %     % read the first image to retrieve its resolution
                        %     I = imread(fileNames(1));
                        
                        % Creation of the collage with the other images
                        for indexmontage = 1 : nummont
                            
                            if indexmontage == nummont
                                
                                clear I
                                
                                % number of images not available to fill the 3*2 collage
                                % imadiff = (nummont * 6) - size(dirOutput,1);
                                
                                for indexfig = (6 * (indexmontage - 1) + 1) : size(dirOutput, 1)
                                    if ispc
                                        I{indexfig - (6 * (indexmontage - 1))} = imread([dirOutput(indexfig).folder '\' dirOutput(indexfig).name]);
                                    else
                                        I{indexfig - (6 * (indexmontage - 1))}  = imread([dirOutput(indexfig).folder '/' dirOutput(indexfig).name]);
                                    end
                                end
                                
                                
                                % Creating fake images to fill the 3x2 image
                                for indexfig = (size(dirOutput,1) + 1) : (6 * indexmontage)
                                    I{indexfig - (6 * (indexmontage - 1))} = 255 * ones(size(I{1, 1}), 'uint8');
                                end
                                
                                Im = [I{1, 1} I{1, 2};
                                    I{1, 3} I{1, 4};
                                    I{1, 5} I{1, 6}];
                                imshow(Im, [], 'InitialMagnification','fit');
                                
                                if ispc
                                    % saveas(f, [handles.FIELDRT.env.userDir '\FIELDRT_Reports\Images\' module '\' struct_name_case{index1} '\' struct_name_temp '_montage_' int2str(indexmontage)] , 'tiff')
                                    % print(f, [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage)] , '-dtiffn', '-r500')
                                    imwrite(Im,[handles.FIELDRT.env.userDir '\FIELDRT_Reports\Images\' module '\' struct_name_case{index1} '\' struct_name_temp '_montage_' int2str(indexmontage) '.tiff'],'tiff')
                                else
                                    % saveas(f, [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage)] , 'tiff')
                                    % print(f, [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage)] , '-dtiffn', '-r500')
                                    imwrite(Im,[handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage) '.tiff'],'tiff')
                                end
                                
                                %             montage(fileNames, 'ThumbnailSize', [(size(I, 1) * 2) (size(I, 2) * 3)], 'Indices', (6 * (indexmontage - 1) + 1) : size(dirOutput,1), 'BackgroundColor', 'white', 'Size', [2 3]);
                                %
                                %             % Saving the montage
                                %             if ispc
                                %                 % saveas(f, [handles.FIELDRT.env.userDir '\FIELDRT_Reports\Images\' module '\' struct_name_case{index1} '\' struct_name_temp '_montage_' int2str(indexmontage)] , 'tiff')
                                %                 print(f, [handles.FIELDRT.env.userDir '\FIELDRT_Reports\Images\' module '\' struct_name_case{index1} '\' struct_name_temp '_montage_' int2str(indexmontage)] , '-dtiffn', '-r500')
                                %             else
                                %                 % saveas(f, [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage)] , 'tiff')
                                %                 print(f, [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage)] , '-dtiffn', '-r500')
                                %             end
                            else
                                % Only for the first structure a collage of four images is
                                % created. This is needed because of the layout of the first page of
                                % each functionality (i.e over/under, OARs...)
                                
                                
                                clear I
                                
                                for indexfig = (6 * (indexmontage - 1) + 1) : (6 * indexmontage)
                                    if ispc
                                        I{indexfig - (6 * (indexmontage - 1))} = imread([dirOutput(indexfig).folder '\' dirOutput(indexfig).name]);
                                    else
                                        I{indexfig - (6 * (indexmontage - 1))}  = imread([dirOutput(indexfig).folder '/' dirOutput(indexfig).name]);
                                    end
                                end
                                
                                Im = [I{1, 1} I{1, 2};
                                    I{1, 3} I{1, 4};
                                    I{1, 5} I{1, 6}];
                                imshow(Im, [], 'InitialMagnification','fit');
                                
                                if ispc
                                    % saveas(f, [handles.FIELDRT.env.userDir '\FIELDRT_Reports\Images\' module '\' struct_name_case{index1} '\' struct_name_temp '_montage_' int2str(indexmontage)] , 'tiff')
                                    % print(f, [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage)] , '-dtiffn', '-r500')
                                    imwrite(Im,[handles.FIELDRT.env.userDir '\FIELDRT_Reports\Images\' module '\' struct_name_case{index1} '\' struct_name_temp '_montage_' int2str(indexmontage) '.tiff'],'tiff')
                                else
                                    % saveas(f, [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage)] , 'tiff')
                                    % print(f, [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage)] , '-dtiffn', '-r500')
                                    imwrite(Im,[handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage) '.tiff'],'tiff')
                                end
                                
                                
                                %             montage(fileNames, 'ThumbnailSize', [(size(I, 1) * 2) (size(I, 2) * 3)], 'Indices', (6 * (indexmontage - 1) + 1) : (6 * indexmontage), 'BackgroundColor', 'white', 'Size', [2 3]);
                                %
                                %             if ispc
                                %                 % saveas(f, [handles.FIELDRT.env.userDir '\FIELDRT_Reports\Images\' module '\' struct_name_case{index1} '\' struct_name_temp '_montage_' int2str(indexmontage)] , 'tiff')
                                %                 print(f, [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage)] , '-dtiffn', '-r500')
                                %             else
                                %                 % saveas(f, [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage)] , 'tiff')
                                %                 print(f, [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage)] , '-dtiffn', '-r500')
                                %             end
                            end
                        end
                    end
                    
                    % Number of figures for each structure
                    indeximagestotal (1, index1) = nummont;
                    
                    % Closing the figure
                    clear f
                end
            end
    end
waitingbar_index = waitingbar_index + 1;
end


% % closing the figure
% clear f

% %% Creating the montage with the created figures
%
% %     % figure for the montage
% %     f = figure('visible','off', 'Position', [screen(1) screen(2) (screen(3) + 500) (screen(4) + 200) ]);
% %
% %     % To have the figure with Minimal White Space
% %     ax = gca;
% %     outerpos = ax.OuterPosition;
% %     ti = ax.TightInset;
% %     left = outerpos(1) + ti(1);
% %     bottom = outerpos(2) + ti(2);
% %     ax_width = outerpos(3) - ti(1) - ti(3);
% %     ax_height = outerpos(4) - ti(2) - ti(4);
% %     ax.Position = [left bottom ax_width ax_height];
%
%     if index1 == 1
%         if ispc
%             fileFolder = [handles.FIELDRT.env.userDir '\FIELDRT_Reports\Images\' module '\' struct_name_case{index1} '\'];
%         else
%             fileFolder = [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/'];
%         end
%
%         dirOutput = dir(fullfile(fileFolder,[struct_name_case{index1} '*.tif']));
%
%         % Only for the first structure a collage of four images is
%         % created. This is needed because of the layout of the first page of
%         % each functionality (i.e over/under, OARs...)
%         dirOutput_temp = dirOutput(1 : 4);
%         fileNames_temp = string({dirOutput_temp.name}); %retrieving filenames
%
%         dirOutput = dirOutput(5 : end);
%         % fileNames = string({dirOutput.name}); %retrieving filenames
%
%
%         % Creation of the first collage with four images for the first
%         % structure
%         clear I
%
%         if size(dirOutput_temp, 1) > 0
%             for indexfig = 1 : size(dirOutput_temp, 1)
%                 if ispc
%                     I{indexfig} = imread([dirOutput_temp(indexfig).folder '\' dirOutput_temp(indexfig).name]);
%                 else
%                     I{indexfig} = imread([dirOutput_temp(indexfig).folder '/' dirOutput_temp(indexfig).name]);
%                 end
%             end
%
%             Im = [I{1, 1} I{1, 2};
%                 I{1, 3} I{1, 4}];
%             % imshow(Im, [], 'InitialMagnification','fit');
%
%             if ispc
%                 % saveas(f, [handles.FIELDRT.env.userDir '\FIELDRT_Reports\Images\' module '\' struct_name_case{index1} '\' struct_name_temp '_montage_' int2str(indexmontage)] , 'tiff')
%                 % print(f, [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage)] , '-dtiffn', '-r500')
%                 imwrite(Im,[handles.FIELDRT.env.userDir '\FIELDRT_Reports\Images\' module '\' struct_name_case{index1} '\' struct_name_temp '_montage_' int2str(1) '.tiff'],'tiff')
%             else
%                 % saveas(f, [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage)] , 'tiff')
%                 % print(f, [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage)] , '-dtiffn', '-r500')
%                 imwrite(Im,[handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(1) '.tiff'],'tiff')
%             end
%         end
%
%
%         % creation of the montage with 6 images
%         if size(dirOutput,1) > 0
%
%             nummont = ceil(size(dirOutput,1) / 6);
%
%             for indexmontage = 1 : nummont
%
%                 if indexmontage == nummont
%
%                     clear I
%
%                     % number of images not available to fill the 3*2 collage
%                     % imadiff = (nummont * 6) - size(dirOutput,1);
%
%                     for indexfig = (6 * (indexmontage - 1) + 1) : size(dirOutput, 1)
%                         if ispc
%                             I{indexfig - (6 * (indexmontage - 1))} = imread([dirOutput(indexfig).folder '\' dirOutput(indexfig).name]);
%                         else
%                             I{indexfig - (6 * (indexmontage - 1))} = imread([dirOutput(indexfig).folder '/' dirOutput(indexfig).name]);
%                         end
%                     end
%
%
%
%
%                     % Creating fake images to fill the 3x2 image
%                     for indexfig = (size(dirOutput,1) + 1) : (6 * indexmontage)
%                         I{indexfig - (6 * (indexmontage - 1))} = 255 * ones(size(I{1, 1}), 'uint8');
%                     end
%
%                     Im = [I{1, 1} I{1, 2};
%                         I{1, 3} I{1, 4};
%                         I{1, 5} I{1, 6}];
%                     imshow(Im, [], 'InitialMagnification','fit');
%
%                     if ispc
%                         % saveas(f, [handles.FIELDRT.env.userDir '\FIELDRT_Reports\Images\' module '\' struct_name_case{index1} '\' struct_name_temp '_montage_' int2str(indexmontage)] , 'tiff')
%                         % print(f, [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage)] , '-dtiffn', '-r500')
%                         imwrite(Im,[handles.FIELDRT.env.userDir '\FIELDRT_Reports\Images\' module '\' struct_name_case{index1} '\' struct_name_temp '_montage_' int2str(indexmontage + 1) '.tiff'],'tiff')
%                     else
%                         % saveas(f, [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage)] , 'tiff')
%                         % print(f, [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage)] , '-dtiffn', '-r500')
%                         imwrite(Im,[handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage + 1) '.tiff'],'tiff')
%                     end
%
%                     %             montage(fileNames, 'ThumbnailSize', [(size(I, 1) * 2) (size(I, 2) * 3)], 'Indices', (6 * (indexmontage - 1) + 1) : size(dirOutput,1), 'BackgroundColor', 'white', 'Size', [2 3]);
%                     %
%                     %             % Saving the montage
%                     %             if ispc
%                     %                 % saveas(f, [handles.FIELDRT.env.userDir '\FIELDRT_Reports\Images\' module '\' struct_name_case{index1} '\' struct_name_temp '_montage_' int2str(indexmontage)] , 'tiff')
%                     %                 print(f, [handles.FIELDRT.env.userDir '\FIELDRT_Reports\Images\' module '\' struct_name_case{index1} '\' struct_name_temp '_montage_' int2str(indexmontage)] , '-dtiffn', '-r500')
%                     %             else
%                     %                 % saveas(f, [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage)] , 'tiff')
%                     %                 print(f, [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage)] , '-dtiffn', '-r500')
%                     %             end
%                 else
%                     % Only for the first structure a collage of four images is
%                     % created. This is needed because of the layout of the first page of
%                     % each functionality (i.e over/under, OARs...)
%
%
%                     clear I
%
%                     for indexfig = (6 * (indexmontage - 1) + 1) : (6 * indexmontage)
%                         if ispc
%                             I{indexfig - (6 * (indexmontage - 1))} = imread([dirOutput(indexfig).folder '\' dirOutput(indexfig).name]);
%                         else
%                             I{indexfig - (6 * (indexmontage - 1))} = imread([dirOutput(indexfig).folder '/' dirOutput(indexfig).name]);
%                         end
%                     end
%
%                     Im = [I{1, 1} I{1, 2};
%                         I{1, 3} I{1, 4};
%                         I{1, 5} I{1, 6}];
%                     imshow(Im, [], 'InitialMagnification','fit');
%
%                     if ispc
%                         % saveas(f, [handles.FIELDRT.env.userDir '\FIELDRT_Reports\Images\' module '\' struct_name_case{index1} '\' struct_name_temp '_montage_' int2str(indexmontage)] , 'tiff')
%                         % print(f, [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage)] , '-dtiffn', '-r500')
%                         imwrite(Im,[handles.FIELDRT.env.userDir '\FIELDRT_Reports\Images\' module '\' struct_name_case{index1} '\' struct_name_temp '_montage_' int2str(indexmontage + 1) '.tiff'],'tiff')
%                     else
%                         % saveas(f, [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage)] , 'tiff')
%                         % print(f, [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage)] , '-dtiffn', '-r500')
%                         imwrite(Im,[handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage + 1) '.tiff'],'tiff')
%                     end
%
%
%                     %             montage(fileNames, 'ThumbnailSize', [(size(I, 1) * 2) (size(I, 2) * 3)], 'Indices', (6 * (indexmontage - 1) + 1) : (6 * indexmontage), 'BackgroundColor', 'white', 'Size', [2 3]);
%                     %
%                     %             if ispc
%                     %                 % saveas(f, [handles.FIELDRT.env.userDir '\FIELDRT_Reports\Images\' module '\' struct_name_case{index1} '\' struct_name_temp '_montage_' int2str(indexmontage)] , 'tiff')
%                     %                 print(f, [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage)] , '-dtiffn', '-r500')
%                     %             else
%                     %                 % saveas(f, [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage)] , 'tiff')
%                     %                 print(f, [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage)] , '-dtiffn', '-r500')
%                     %             end
%                 end
%             end
%         end
%
%         indeximagestotal (index1, 1) = nummont + 1;
%
%     else
%         % creation of the collage for the other structures
%         if ispc
%             fileFolder = [handles.FIELDRT.env.userDir '\FIELDRT_Reports\Images\' module '\' struct_name_case{index1} '\'];
%         else
%             fileFolder = [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/'];
%         end
%
%         dirOutput = dir(fullfile(fileFolder,[struct_name_case{index1} '*.tif']));
%
%         dirOutput = dirOutput(5 : end);
%         fileNames = string({dirOutput.name}); %retrieving filenames
%
%         if size(dirOutput,1) > 0
%             % Numbers of montages --> we want 6 images per page
%             nummont = ceil(size(dirOutput,1) / 6);
%
%             %     % read the first image to retrieve its resolution
%             %     I = imread(fileNames(1));
%
%             % Creation of the collage with the other images
%             for indexmontage = 1 : nummont
%
%                 if indexmontage == nummont
%
%                     clear I
%
%                     % number of images not available to fill the 3*2 collage
%                     % imadiff = (nummont * 6) - size(dirOutput,1);
%
%                     for indexfig = (6 * (indexmontage - 1) + 1) : size(dirOutput, 1)
%                         if ispc
%                             I{indexfig - (6 * (indexmontage - 1))} = imread([dirOutput(indexfig).folder '\' dirOutput(indexfig).name]);
%                         else
%                             I{indexfig - (6 * (indexmontage - 1))}  = imread([dirOutput(indexfig).folder '/' dirOutput(indexfig).name]);
%                         end
%                     end
%
%
%                     % Creating fake images to fill the 3x2 image
%                     for indexfig = (size(dirOutput,1) + 1) : (6 * indexmontage)
%                         I{indexfig - (6 * (indexmontage - 1))} = 255 * ones(size(I{1, 1}), 'uint8');
%                     end
%
%                     Im = [I{1, 1} I{1, 2};
%                         I{1, 3} I{1, 4};
%                         I{1, 5} I{1, 6}];
%                     imshow(Im, [], 'InitialMagnification','fit');
%
%                     if ispc
%                         % saveas(f, [handles.FIELDRT.env.userDir '\FIELDRT_Reports\Images\' module '\' struct_name_case{index1} '\' struct_name_temp '_montage_' int2str(indexmontage)] , 'tiff')
%                         % print(f, [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage)] , '-dtiffn', '-r500')
%                         imwrite(Im,[handles.FIELDRT.env.userDir '\FIELDRT_Reports\Images\' module '\' struct_name_case{index1} '\' struct_name_temp '_montage_' int2str(indexmontage) '.tiff'],'tiff')
%                     else
%                         % saveas(f, [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage)] , 'tiff')
%                         % print(f, [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage)] , '-dtiffn', '-r500')
%                         imwrite(Im,[handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage) '.tiff'],'tiff')
%                     end
%
%                     %             montage(fileNames, 'ThumbnailSize', [(size(I, 1) * 2) (size(I, 2) * 3)], 'Indices', (6 * (indexmontage - 1) + 1) : size(dirOutput,1), 'BackgroundColor', 'white', 'Size', [2 3]);
%                     %
%                     %             % Saving the montage
%                     %             if ispc
%                     %                 % saveas(f, [handles.FIELDRT.env.userDir '\FIELDRT_Reports\Images\' module '\' struct_name_case{index1} '\' struct_name_temp '_montage_' int2str(indexmontage)] , 'tiff')
%                     %                 print(f, [handles.FIELDRT.env.userDir '\FIELDRT_Reports\Images\' module '\' struct_name_case{index1} '\' struct_name_temp '_montage_' int2str(indexmontage)] , '-dtiffn', '-r500')
%                     %             else
%                     %                 % saveas(f, [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage)] , 'tiff')
%                     %                 print(f, [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage)] , '-dtiffn', '-r500')
%                     %             end
%                 else
%                     % Only for the first structure a collage of four images is
%                     % created. This is needed because of the layout of the first page of
%                     % each functionality (i.e over/under, OARs...)
%
%
%                     clear I
%
%                     for indexfig = (6 * (indexmontage - 1) + 1) : (6 * indexmontage)
%                         if ispc
%                             I{indexfig - (6 * (indexmontage - 1))} = imread([dirOutput(indexfig).folder '\' dirOutput(indexfig).name]);
%                         else
%                             I{indexfig - (6 * (indexmontage - 1))}  = imread([dirOutput(indexfig).folder '/' dirOutput(indexfig).name]);
%                         end
%                     end
%
%                     Im = [I{1, 1} I{1, 2};
%                         I{1, 3} I{1, 4};
%                         I{1, 5} I{1, 6}];
%                     imshow(Im, [], 'InitialMagnification','fit');
%
%                     if ispc
%                         % saveas(f, [handles.FIELDRT.env.userDir '\FIELDRT_Reports\Images\' module '\' struct_name_case{index1} '\' struct_name_temp '_montage_' int2str(indexmontage)] , 'tiff')
%                         % print(f, [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage)] , '-dtiffn', '-r500')
%                         imwrite(Im,[handles.FIELDRT.env.userDir '\FIELDRT_Reports\Images\' module '\' struct_name_case{index1} '\' struct_name_temp '_montage_' int2str(indexmontage) '.tiff'],'tiff')
%                     else
%                         % saveas(f, [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage)] , 'tiff')
%                         % print(f, [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage)] , '-dtiffn', '-r500')
%                         imwrite(Im,[handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage) '.tiff'],'tiff')
%                     end
%
%
%                     %             montage(fileNames, 'ThumbnailSize', [(size(I, 1) * 2) (size(I, 2) * 3)], 'Indices', (6 * (indexmontage - 1) + 1) : (6 * indexmontage), 'BackgroundColor', 'white', 'Size', [2 3]);
%                     %
%                     %             if ispc
%                     %                 % saveas(f, [handles.FIELDRT.env.userDir '\FIELDRT_Reports\Images\' module '\' struct_name_case{index1} '\' struct_name_temp '_montage_' int2str(indexmontage)] , 'tiff')
%                     %                 print(f, [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage)] , '-dtiffn', '-r500')
%                     %             else
%                     %                 % saveas(f, [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage)] , 'tiff')
%                     %                 print(f, [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp '_montage_' int2str(indexmontage)] , '-dtiffn', '-r500')
%                     %             end
%                 end
%             end
%         end
%
%         % Number of figures for each structure
%         indeximagestotal (index1, 1) = nummont;
%
%         % Closing the figure
%         clear f
%     end
%