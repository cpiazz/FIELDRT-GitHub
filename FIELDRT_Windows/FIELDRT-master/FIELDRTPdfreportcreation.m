function FIELDRTPdfreportcreation(varargin)
% Creation of the pdf report (with only selected images)

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

% indexcase = handles.Case; % Case selected


% Waiting bar
wb = waitbar(0, 'Creating the images for the report ...');
waitingbar_index = 1;
% if indexcase == 11 % Oesophagus Case1
%     waitingbar_struct_name_case = 13; % 5 structs for Over/under, 2 structs for OARs, 5 structs for min and max, 1 fake structure to avoid to have 
%     % the waiting bar filled during the creation of the images for the last structure for min and max    
% end
% 
% 
% if indexcase == 31 % Prostate Case1,
%     waitingbar_struct_name_case = 7; % 2 structs for Over/under, 2 structs for OARs, 2 structs for min and max, 1 fake structure to avoid to have 
%     % the waiting bar filled during the creation of the images for the last structure for min and max    
% end

% the total lenght of the waiting bar is defined by:
% - number of structures for the over/under module 
% - number of structures for the OARs module 
% - number of structures for the min and max module 
% - 1 fake structure to avoid to have the waiting bar filled during the creation of the images for the last structure for min and max    
waitingbar_struct_name_case = (size(handles.FIELDRTGSCases(handles.Case(1,1)).StructOverunder(handles.Case(1,2)).Cases, 2) + ...
    size(handles.FIELDRTGSCases(handles.Case(1,1)).StructOARs(handles.Case(1,2)).Cases, 2) + ...
    size(handles.FIELDRTGSCases(handles.Case(1,1)).StructMinmax(handles.Case(1,2)).Cases, 2) + ...
    1);

%% Creating the images for each module
for indexmodule = 1 : 3
    % case1 = over/under contoured regions
    % case2 = OARs
    % case3 = max/min
    
    switch indexmodule
        
        case 1 % Over and under module   
            
%             if indexcase == 11 % Oesophagus Case1
%                 struct_name_case = {'GTV', 'CTVA', 'CTVB', 'CTVC', 'PTV'};
%                 struct_name_case_OARs = {}; % No OARs needed for this module            
%             end
%             
%             if indexcase == 31 % Prostate Case1
%                 struct_name_case = {'CTVp', 'CTVpsv'};
%                 struct_name_case_OARs = {}; % No OARs needed for this module            
%             end

            struct_name_case = handles.FIELDRTGSCases(handles.Case(1,1)).StructOverunder(handles.Case(1,2)).Cases;
            struct_name_case_OARs = {}; % No OARs needed for this module  
            
            module = 'Overunder';
            
            
            %% Finding the scan range for each structure
            % The range is determined by the lowest/highest scan in which
            % the GS or the user contour are appearing
            for index1 = 1 : size(struct_name_case, 2)
                % Creating a dummy folder in the report folder to save there all the images for the report
                if ispc
                    mkdir([handles.FIELDRT.env.userDir '\FIELDRT_Reports\Images\' module '\' struct_name_case{index1} '\']);
                else
                    mkdir([handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/']);
                end
                
                range = [];
                struct_name_temp = struct_name_case{index1};
                % Retrieving information about the selected structure
                structnum  = getStructNum(struct_name_temp, handles.planC, indexS);
                
                if iscolumn(handles.planC{1, 4}(structnum).contour)
                    for index2 = 1 : size(handles.planC{1, 4}(structnum).contour, 1)
                        if size(handles.planC{1, 4}(structnum).contour(index2).segments, 2) ~= 0
                            pointsM = handles.planC{handles.structfieldnum}(structnum).contour(index2).segments.points;
                            if isempty(pointsM) % nothing happens
                                
                            else
                                range = [range; index2];
                            end
                        end
                    end
                else
                    for index2 = 1 : size(handles.planC{1, 4}(structnum).contour, 2)
                        if size(handles.planC{1, 4}(structnum).contour(index2).segments, 2) ~= 0
                            pointsM = handles.planC{handles.structfieldnum}(structnum).contour(index2).segments.points;
                            if isempty(pointsM) % nothing happens
                                
                            else
                                range = [range; index2];
                            end
                        end
                    end
                end
                
                
                range_GS = [];
                struct_name_temp_GS = [struct_name_case{index1} '_GS'];
                % Retrieving information about the selected structure
                structnum_GS  = getStructNum(struct_name_temp_GS, handles.planC, indexS);
               
                if iscolumn(handles.planC{1, 4}(structnum_GS).contour)
                    for index2 = 1 : size(handles.planC{1, 4}(structnum_GS).contour, 1)
                        if size(handles.planC{1, 4}(structnum_GS).contour(index2).segments, 2) ~= 0
                            pointsM = handles.planC{handles.structfieldnum}(structnum_GS).contour(index2).segments.points;
                            if isempty(pointsM) % nothing happens
                                
                            else
                                range_GS = [range_GS; index2];
                            end
                        end
                    end
                else
                    for index2 = 1 : size(handles.planC{1, 4}(structnum_GS).contour, 2)
                        if size(handles.planC{1, 4}(structnum_GS).contour(index2).segments, 2) ~= 0
                            pointsM = handles.planC{handles.structfieldnum}(structnum_GS).contour(index2).segments.points;
                            if isempty(pointsM) % nothing happens
                                
                            else
                                range_GS = [range_GS; index2];
                            end
                        end
                    end
                end
                
                % Computing the lowest slice
                if min(range) < min(range_GS)
                    range_total_min = min(range);
                else
                    range_total_min = min(range_GS);
                end
                
                if max(range) > max(range_GS)
                    range_total_max = max(range);
                else
                    range_total_max = max(range_GS);
                end
                
                % Final range
                range_total(index1, :) = [range_total_min range_total_max];
            end
            
            %% Creating the traffic light conformity image and flagging the
            % images to show (the ones with 2D Jaccard < 0.50)
            if ispc
                mkdir([handles.FIELDRT.env.userDir '\FIELDRT_Reports\Images\' module '\Stats\']);
            else
                mkdir([handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/Stats/']);
            end
            
            [Overunder_imagestoshow, handles.jaccard2Doutput] = FIELDRTstatsfiguresaving(varargin, range_total, struct_name_case, module);
            
            %% Creating and saving the figure to save in the pdf report
            [indeximagestotalmodule(indexmodule).numberofmontage, waitingbar_index] = FIELDRTPdffiguresaving(varargin, range_total, struct_name_case, module, struct_name_case_OARs, Overunder_imagestoshow, wb, waitingbar_index, waitingbar_struct_name_case);
            
        case 2 % OARs module
            
            %             indexcase = handles.Case; % Case selected
            %
            %             if indexcase == 11 % Oesophagus Case1
            %                 struct_name_case = {'GTV', 'CTVB'};
            %                 struct_name_case_OARs = {'Vertebra_GS'; 'Aorta_GS'; 'Right lung_GS'; 'Pericardium/great vessels_GS'; ...
            %                     'Liver_GS'; 'Stomach_GS'; 'Azygous vein_GS'; 'Left main bronchus_GS'; 'Left lung_GS'};
            %             end
            %
            %             if indexcase == 31 % Prostate Case1
            %                 struct_name_case = {'CTVp', 'CTVpsv'};
            %                 struct_name_case_OARs = {'Bladder'; 'Bowel'; 'Lt FemHead'; 'Penile bulb'; 'Rectum'; 'Rt_FemHead'};
            %             end
            
            struct_name_case = handles.FIELDRTGSCases(handles.Case(1,1)).StructOARs(handles.Case(1,2)).Cases;
            struct_name_case_OARs = handles.FIELDRTGSCases(handles.Case(1,1)).OARs(handles.Case(1,2)).Cases; % No OARs needed for this module
            
            module = 'OARs';
            
            %% Finding the scan range for each structure
            % The range is determined by the lowest/highest scan in which
            % the GS or the user contour are appearing
            for index1 = 1 : size(struct_name_case, 2)
                % Creating a dummy folder in the report folder to save there all the images for the report
                if ispc
                    mkdir([handles.FIELDRT.env.userDir '\FIELDRT_Reports\Images\' module '\' struct_name_case{index1} '\']);
                else
                    mkdir([handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/']);
                end
                
                range = [];
                struct_name_temp = struct_name_case{index1};
                % Retrieving information about the selected structure
                structnum  = getStructNum(struct_name_temp, handles.planC, indexS);
                if iscolumn(handles.planC{1, 4}(structnum).contour)
                    for index2 = 1 : size(handles.planC{1, 4}(structnum).contour, 1)
                        if size(handles.planC{1, 4}(structnum).contour(index2).segments, 2) ~= 0
                            pointsM = handles.planC{handles.structfieldnum}(structnum).contour(index2).segments.points;
                            if isempty(pointsM) % nothing happens
                                
                            else
                                range = [range; index2];
                            end
                        end
                    end
                else
                    for index2 = 1 : size(handles.planC{1, 4}(structnum).contour, 2)
                        if size(handles.planC{1, 4}(structnum).contour(index2).segments, 2) ~= 0
                            pointsM = handles.planC{handles.structfieldnum}(structnum).contour(index2).segments.points;
                            if isempty(pointsM) % nothing happens
                                
                            else
                                range = [range; index2];
                            end
                        end
                    end
                end
                
                
                range_GS = [];
                struct_name_temp_GS = [struct_name_case{index1} '_GS'];
                % Retrieving information about the selected structure
                structnum_GS  = getStructNum(struct_name_temp_GS, handles.planC, indexS);
                if iscolumn(handles.planC{1, 4}(structnum_GS).contour)
                    for index2 = 1 : size(handles.planC{1, 4}(structnum_GS).contour, 1)
                        if size(handles.planC{1, 4}(structnum_GS).contour(index2).segments, 2)  ~= 0
                            pointsM = handles.planC{handles.structfieldnum}(structnum_GS).contour(index2).segments.points;
                            if isempty(pointsM) % nothing happens
                                
                            else
                                range_GS = [range_GS; index2];
                            end
                        end
                    end
                else
                    for index2 = 1 : size(handles.planC{1, 4}(structnum_GS).contour, 2)
                        if size(handles.planC{1, 4}(structnum_GS).contour(index2).segments, 2)  ~= 0
                            pointsM = handles.planC{handles.structfieldnum}(structnum_GS).contour(index2).segments.points;
                            if isempty(pointsM) % nothing happens
                                
                            else
                                range_GS = [range_GS; index2];
                            end
                        end
                    end
                end
                
                % Computing the lowest slice
                if min(range) < min(range_GS)
                    range_total_min = min(range);
                else
                    range_total_min = min(range_GS);
                end
                
                if max(range) > max(range_GS)
                    range_total_max = max(range);
                else
                    range_total_max = max(range_GS);
                end
                
                % Final range
                range_total(index1, :) = [range_total_min range_total_max];
            end
            
            %% Identifying the slices in which the US contour is going too much into the OARs
            % Creating a matrix so to flag the slices to show
            OARs_imagestoshow = uint8(zeros([size(struct_name_case, 2) size(cell2mat(handles.jaccard2Doutput), 2)]));
            
            % For each OAR and for each structure, checking if the slice has to be
            % flagged
            for index1_OARs = 1 : size(OARs_imagestoshow, 1) % structure iterator
                for index2_OARs = range_total(index1_OARs, 1) : range_total(index1_OARs, 2)  % slice iterator
                    indexout = 0; % just a flag to stop the following
                    while OARs_imagestoshow(index1_OARs, index2_OARs) == 0 && indexout == 0
                        for index3_OARs = 1 : size(handles.OARsdummy_finalmask, 2) % OARs iterator
                            tempmask = handles.OARsdummy_finalmask{index1_OARs, index3_OARs}(:, :, index2_OARs); % Estracting the slice
                            clear CC
                            CC = bwconncomp(tempmask); % connectivity = 8 (default)
                            if CC.NumObjects > 0
                                % TT = regionprops(CC, 'Centroid', 'Area');
                                % for indexNumObjects = 1 : CC.NumObjects
                                % if TT(indexNumObjects, 1).Area > 30 % We don't want to annotate just a single point
                                % % flag current slice
                                % OARs_imagestoshow(index1_OARs, index2_OARs) = 1;
                                % end
                                % end
                                OARs_imagestoshow(index1_OARs, index2_OARs) = 1;
                            end
                        end
                        indexout = 1;
                    end
                end
            end
            
             %%  Creating and saving the figure to save in the pdf report
            [indeximagestotalmodule(indexmodule).numberofmontage, waitingbar_index] = FIELDRTPdffiguresaving(varargin, range_total, struct_name_case, module, struct_name_case_OARs, OARs_imagestoshow, wb, waitingbar_index, waitingbar_struct_name_case);
            
        case 3 % Min and max module
                      
            %             indexcase = handles.Case; % Case selected
            %
            %             if indexcase == 11 % Oesophagus Case1
            %                 struct_name_case = {'GTV', 'CTVA', 'CTVB', 'CTVC', 'PTV'};
            %                 struct_name_case_OARs = {}; % No OARs needed for this module
            %             end
            %
            %             if indexcase == 31 % Prostate Case1
            %                 struct_name_case = {'CTVp', 'CTVpsv'};
            %                 struct_name_case_OARs = {}; % No OARs needed for this module
            %             end
            struct_name_case = handles.FIELDRTGSCases(handles.Case(1,1)).StructMinmax(handles.Case(1,2)).Cases;
            struct_name_case_OARs = {}; % No OARs needed for this module
                        
            module = 'Minmax';
            
            %% Finding the scan range for each structure
            % The range is determined by the lowest/highest scan in which
            % the GS or the user contour are appearing
            for index1 = 1 : size(struct_name_case, 2)
                % Creating a dummy folder in the report folder to save there all the images for the report
                if ispc
                    mkdir([handles.FIELDRT.env.userDir '\FIELDRT_Reports\Images\' module '\' struct_name_case{index1} '\']);
                else
                    mkdir([handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/']);
                end
                
                range = [];
                struct_name_temp = struct_name_case{index1};
                % Retrieving information about the selected structure
                structnum  = getStructNum(struct_name_temp, handles.planC, indexS);
                if iscolumn(handles.planC{1, 4}(structnum).contour)
                    for index2 = 1 : size(handles.planC{1, 4}(structnum).contour, 1)
                        if size(handles.planC{1, 4}(structnum).contour(index2).segments, 2) ~= 0
                            pointsM = handles.planC{handles.structfieldnum}(structnum).contour(index2).segments.points;
                            if isempty(pointsM) % nothing happens
                                
                            else
                                range = [range; index2];
                            end
                        end
                    end
                else
                    for index2 = 1 : size(handles.planC{1, 4}(structnum).contour, 2)
                        if size(handles.planC{1, 4}(structnum).contour(index2).segments, 2) ~= 0
                            pointsM = handles.planC{handles.structfieldnum}(structnum).contour(index2).segments.points;
                            if isempty(pointsM) % nothing happens
                                
                            else
                                range = [range; index2];
                            end
                        end
                    end
                end
                
                
                range_GS = [];
                struct_name_temp_GS = [struct_name_case{index1} '_GS'];
                % Retrieving information about the selected structure
                structnum_GS  = getStructNum(struct_name_temp_GS, handles.planC, indexS);
                if iscolumn(handles.planC{1, 4}(structnum_GS).contour)
                    for index2 = 1 : size(handles.planC{1, 4}(structnum_GS).contour, 1)
                        if size(handles.planC{1, 4}(structnum_GS).contour(index2).segments, 2)  ~= 0
                            pointsM = handles.planC{handles.structfieldnum}(structnum_GS).contour(index2).segments.points;
                            if isempty(pointsM) % nothing happens
                                
                            else
                                range_GS = [range_GS; index2];
                            end
                        end
                    end
                else
                    for index2 = 1 : size(handles.planC{1, 4}(structnum_GS).contour, 2)
                        if size(handles.planC{1, 4}(structnum_GS).contour(index2).segments, 2)  ~= 0
                            pointsM = handles.planC{handles.structfieldnum}(structnum_GS).contour(index2).segments.points;
                            if isempty(pointsM) % nothing happens
                                
                            else
                                range_GS = [range_GS; index2];
                            end
                        end
                    end
                end
                
                % Computing the lowest slice
                if min(range) < min(range_GS)
                    range_total_min = min(range);
                else
                    range_total_min = min(range_GS);
                end
                
                if max(range) > max(range_GS)
                    range_total_max = max(range);
                else
                    range_total_max = max(range_GS);
                end
                
                % Final range
                range_total(index1, :) = [range_total_min range_total_max];
            end
            
            %% Identifying the slices to show (currently the min and max contour are not showed in the PDF report)
            % Creating a matrix so to flag the slices to show
            Minmax_imagestoshow = uint8(zeros(size(cell2mat(handles.jaccard2Doutput))));          
            
            %% Creating and saving the figure to save in the pdf report
            [indeximagestotalmodule(indexmodule).numberofmontage, waitingbar_index] = FIELDRTPdffiguresaving(varargin, range_total, struct_name_case, module, struct_name_case_OARs, Minmax_imagestoshow, wb, waitingbar_index, waitingbar_struct_name_case);
            
    end
end

% Updating the waiting bar
wb = waitbar(((waitingbar_index + 1) / waitingbar_struct_name_case), wb);

    
close(wb);

%% Report creation

% Waiting bar
wb = waitbar(0, 'Creating the pdf report ...');

if isdeployed == 1 % checking whether the MATLAB or deployed version is running
    % 0 --> MATLAB version
    % 1 --> deployed version
                            
    % To enable someone who does not have MATLAB installed to run the
    % report generation part
    makeDOMCompilable(); 
end

% Import report API classes (optional)
import mlreportgen.report.*
import mlreportgen.dom.*

% Removing potential space in the name of the attempt
modifiedStr = strrep(handles.AttemptInformation.FileName(1, 1 : (end - 4)), ' ',  '_');
% Add report container (required)
% Creating the report folder
if ispc
    rpt = Report([handles.FIELDRT.env.userDir  '\FIELDRT_Reports\' modifiedStr],'pdf');
else
    rpt = Report([handles.FIELDRT.env.userDir  '/FIELDRT_Reports/' modifiedStr],'pdf');
end

currentLayout.PageMargins.Top = '0.25in';
currentLayout.PageMargins.Header = '0.15in';

% Add content to container (required)
% Types of content added here: title
% page and table of contents reporters
% Adding an empty section with t
              
titlepg = TitlePage;

if ispc
    titlepg.Image = which([handles.FIELDRT.env.userDir '\FIELDRT-master\Utilities\FIELDRT_logo.png']);    
    imgtitle = Image(which([handles.FIELDRT.env.userDir '\FIELDRT-master\Utilities\FIELDRT_logo.png']));
else
    % titlepg.Image = which([handles.FIELDRT.env.userDir '/FIELDRT-master/Utilities/FIELDRT_logo.png']);
    imgtitle = Image(which([handles.FIELDRT.env.userDir '/FIELDRT-master/Utilities/FIELDRT_logo.png']));
end

imgtitle.Style = {ScaleToFit(false), HAlign('center'), Height('3in')};
% img.Style = {HAlign('center')}; %, Height('13in')};
                
titlepg.Image = imgtitle;

supertitle1 = Text('FIELD');

supertitle2 = Text('RT');
supertitle2.Style = {VerticalAlign('superscript'),FontSize('17pt')};

supertitle3 = Text(' report');

titlepg.Title = [supertitle1 supertitle2 supertitle3];
titlepg.Author = ['"' handles.AttemptInformation.FileName(1, 1 : (end - 4)) '" attempt'];

add(rpt,titlepg);
t = TableOfContents;
add(rpt,TableOfContents);

% Add content to report sections (optional)
% Text and formal image added to chapter
% Chapter 1 - Over/under contoured regions
% ch1 = Chapter; 
% ch1.Title = 'Over/under contoured regions'; 
% sec1 = Section; 
% sec1.Title = 'Over and under regions obtained by the comparison of the user outlined structures with the reference ones'; 
% add(ch1, sec1);
% add(rpt, ch1);


% chTitle = Heading1('Over/under contoured regions');
% add(chap,'Over and under regions obtained by the comparison of the user outlined structures with the reference ones'); 

% Style for each chapter
chTitle = Heading1();
chTitle.Style = {CounterInc('sect1'),...
    WhiteSpace('preserve')...
    Color('black'),...
    Bold, FontSize('20pt'), FontFamily('Calibri')};

% Style for each section (description of the chapter)
sectTitle = Heading3();
sectTitle.Style = {CounterInc('sect2'),...
    WhiteSpace('preserve'), Bold(false), LineSpacing('18pt'), FirstLineIndent('0pt'), ...
    Color('black'), FontSize('12pt'), HAlign('justify'), FontFamily('Calibri')};

% Style for each paragraph (structure)
parTitle = Heading2();
parTitle.Style = {CounterInc('sect3'),...
    WhiteSpace('preserve'), Bold(true), FirstLineIndent('0pt')...
    Color('black'), HAlign('left'), FontSize('16pt'), FontFamily('Calibri'),};

% Style for each structure for the statistics
parTitle_stat = Heading2();
parTitle_stat.Style = { WhiteSpace('preserve'), Bold(true), FirstLineIndent('6pt')...
    Color('black'), HAlign('left'), FontSize('16pt'), FontFamily('Calibri')};

% Style for the title of the structures
parTitle_image = Heading3();
parTitle_image.Style = { WhiteSpace('preserve'), Bold(true), FirstLineIndent('6pt')...
    Color('black'), HAlign('left'), FontSize('16pt'), FontFamily('Calibri')};

% Style for the text in min/max section with FIELDRT
parTitle_minmaxFIELDRT = Heading3();
parTitle_minmaxFIELDRT.Style = {WhiteSpace('preserve'), Bold(false), LineSpacing('18pt'), FirstLineIndent('0pt'), ...
    Color('black'), FontSize('12pt'), HAlign('justify'), FontFamily('Calibri')};

%% Creating the images for each module
for indexmodule = 1 : 3
    % Over and under regions
    
    if indexmodule == 1
        
        % Updating the waiting bar
        wb = waitbar((indexmodule / 8), wb);
        
        module = 'Overunder';
        modulenum = 1;
        
        struct_name_case = [];
        % Retrieving structures for this module
        struct_name_case = handles.FIELDRTGSCases(handles.Case(1,1)).StructOverunder(handles.Case(1,2)).Cases;
        
        %% First page with Module description and what it is going to be showed
        chtitleover = clone(chTitle);
        append(chtitleover, AutoNumber('sect1'));
        append(chtitleover,'. ');
        append(chtitleover,'Over/under contoured regions');
        ch1 = Chapter('Title', chtitleover);
        
        paratext = ['This section shows the over and under contoured regions obtained by the comparison between '...
                      'the user outlined structures and the the reference ones.'];
                  
        paratext1 = 'An under contoured region (showed as a cyan contour) is a region present in the reference contour but not outlined by the user.';
        paratext2 = 'An over contoured region (showed as a red contour) is a region outlined by the user but not present in the reference contour.';
        paratext3 = 'Quantitative evaluation of the user outline was performed in terms of volume, 2D and 3D Jaccard similarity index.';
        paratext4 = ['The Jaccard similarity index (sometimes called the Jaccard similarity coefficient) compares members for two sets to see ' ...
                      'which members are shared and which are distinct. It is a measure of similarity for the two sets of data, with a range from '... 
                      '0% to 100%. The higher the percentage, the more similar the two populations. It is computed as follow:'];
                 
        para1 = Text(paratext1);
        para1.Style = sectTitle.Style;
        
        para2 = Text(paratext2);
        para2.Style = sectTitle.Style;
        
        para3 = Text(paratext3);
        para3.Style = sectTitle.Style;
        
        para4 = Text(paratext4);
        para4.Style = sectTitle.Style;
            
        % Jaccard image
        if ispc
            img = Image([handles.FIELDRT.env.userDir '\FIELDRT-master\Utilities\Jaccard_figure.png']);
        else 
            img = Image([handles.FIELDRT.env.userDir '/FIELDRT-master/Utilities/Jaccard_figure.png']);
        end
        img.Style = {ScaleToFit(false), HAlign('center'), Height('3.5in')};
        % img.Style = {HAlign('center')}; %, Height('13in')};
        
        paratext5 = ['In the following pages a summary of the global and local (only 2D Jaccard) quantitative perfomance of the user will be presented. '...
                     'In addition, slices where the 2D Jaccard was below a threshold (0.50) are showed.'];
                 
        para5 = Text(paratext5);
        para5.Style = sectTitle.Style;
        
        % Line break
        l = LineBreak();
        
        sec0 = Section();
        % Adding everything
        add(sec0, para1);
        add(sec0, para2);
        add(sec0, para3);
        add(sec0, para4);
        add(sec0,l); % Adding a line break
        add(sec0, img);
        add(sec0,l); % Adding a line break
        add(sec0, para5);
                
        add(ch1,sec0);
        
%         % Create and append the page break object
%         br = PageBreak();
%         append(ch1,br);
        
        %% Second page with the quantitative evaluation (table with general quantitative evaluation and traffic light
%         structtitleover = clone(parTitle);
%         append(structtitleover, int2str(modulenum));
%         append(structtitleover,'. ');
%         append(structtitleover, '1');
%         append(structtitleover,'. Quantitative evaluation');
%         
%         sec11 = Section('Title', structtitleover);   
        sec12 = Section();  
        
        % Computing the mean 2D Jaccard for each structure and creating the
        % matrix with all the data to use for the table
        for indexjaccard = 1 : size(struct_name_case, 2)
            meanjac2D(indexjaccard, 1) = mean(cell2mat(handles.jaccard2Doutput(indexjaccard, (range_total(indexjaccard, 1) : range_total(indexjaccard, 2)))));
            testData_raw(indexjaccard, 1) = {handles.volumestat_GS{indexjaccard, 1}.vol};
            testData_raw(indexjaccard, 2) = {handles.volumestat_US{indexjaccard, 1}.vol};
            testData_raw(indexjaccard, 3) = {handles.jaccard3Doutput{indexjaccard, 1}};
            % testData_raw(indexjaccard, 4) = {meanjac2D(indexjaccard, 1)};
            
            % An interger copy of the variable testData_raw
            testData_raw_int(indexjaccard, 1) = handles.volumestat_GS{indexjaccard, 1}.vol;
            testData_raw_int(indexjaccard, 2) = handles.volumestat_US{indexjaccard, 1}.vol;
            testData_raw_int(indexjaccard, 3) = handles.jaccard3Doutput{indexjaccard, 1};
            % testData_raw_int(indexjaccard, 4) = meanjac2D(indexjaccard, 1);
        end

        % Obtain cell array size
        [nrows,ncols] = size(testData_raw);
                
        % Preallocate memory for cell array
        testData{nrows,ncols} = [];
        
        % Convert all values to strings to control number of
        % decimal places displayed
        testData = testData_raw;
        idx = cellfun(@isnumeric, testData_raw(:));
        testData(idx) = cellfun(@(x){sprintf('%.2f', x)}, testData_raw(idx));
        
        flag = 0; %flag if the Jaccard is below 0.50
        % Set color of results column text items
        for i = 1 : nrows
            for j = 1 : ncols
                d = string(testData(i,j));
                p = Paragraph(d);
                if j == 3 % || j == 4
                    if testData_raw_int(i,j) < 0.50
                        p.BackgroundColor = 'red';
                        flag = 1;
                    end
                    
                    if testData_raw_int(i,j) < 0.75 && testData_raw_int(i,j) >= 0.50
                        p.BackgroundColor = 'yellow';
                        flag = 1;
                    end
                end
                p.HAlign = 'center';
                testData(i,j) = {p};
            end
        end
        
        testData = [struct_name_case' testData];
        
        % Create and format table
        table = FormalTable({'    ', 'Reference volume [ml]', 'User volume [ml]', '3D Jaccard'}, testData);
        table.RowSep = 'Solid';
        table.ColSep = 'Solid';
        table.Border = 'Solid';
        table.TableEntriesInnerMargin = '8px';
        table.Header.Style = {Bold(), HAlign('center')};
        table.OuterLeftMargin = '12px';
        table.HAlign = 'center';
        table.Width = '180mm';
        
        % table.HAlign = 'right';
        % table.Width = '100px';
        
        structnumb = Text(['Jaccard values below the threshold (0.50) are highlighted in red. Jaccard values between 0.50 and 0.75 are highlighted in yellow']);
        structnumb.Style = {WhiteSpace('preserve'), Bold(false), FirstLineIndent('10pt')...
            Color('black'), HAlign('center'), FontSize('10pt'), FontFamily('Calibri')};
        
        add(sec12,l); % Adding a line break
        add(sec12,table);
        % add(sec13,l); % Adding a line break
        add(sec12, structnumb);        
        
        % Loading the traffic light image
        if ispc% Folder with the images
            fileFolder = [handles.FIELDRT.env.userDir '\FIELDRT_Reports\Images\' module '\Stats\'];
            dirOutput = dir(fullfile(fileFolder, '*_stats.tif'));
        else % Folder with the images
            fileFolder = [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/Stats/'];
            dirOutput = dir(fullfile(fileFolder, '*_stats.tif'));
        end        
        
        % Defining the style for the image
        img = Image([fileFolder dirOutput(1).name]);
        img.Style = {ScaleToFit(false), HAlign('center'), Height('2.1in')};
        % img.Style = {HAlign('center')}; %, Height('13in')};
        
        % Loading the traffic light image legend
        if ispc% Folder with the images
            imglegend = Image([handles.FIELDRT.env.userDir '\FIELDRT-master\Utilities\Traffic_light_legend.png']);
        else
            imglegend = Image([handles.FIELDRT.env.userDir '/FIELDRT-master/Utilities/Traffic_light_legend.png']);
        end        
        
        % Defining the style for the image
        imglegend.Style = {ScaleToFit(false), HAlign('center'), Height('0.3in')};
        % img.Style = {HAlign('center')}; %, Height('13in')};
        
        % footnote
        structnumb1 = Text('Statistics - (1/1)');
        structnumb1.Style = sectTitle.Style;
        structnumb1.Style = {VAlign('bottom'), HAlign('right')};
        
        % Adding everything to the section
        add(sec12, img);
        
        add(sec12, imglegend);
        
        add(sec12, structnumb1);        
        
        structtitleover = clone(parTitle);
        append(structtitleover, int2str(modulenum));
        append(structtitleover,'. ');
        append(structtitleover, '1');
        append(structtitleover,'. Quantitative evaluation');
        
        ch12 = Chapter('Title', structtitleover);
        ch12.Layout.Landscape = true; % Layout of the page containing the chapter
        add(ch12, sec12); % adding the section to the chapter
        
        %% Third page with the flagged images        
        structtitleover = clone(parTitle);
        append(structtitleover, int2str(modulenum));
        append(structtitleover,'. ');
        append(structtitleover, '2');
        append(structtitleover,'. Qualitative evaluation');
        ch13 = Chapter('Title', structtitleover);        
        
        ch13.Layout.Landscape = false; % Layout of the page containing the chapter
                                
        indeximagestotalmodule_int = indeximagestotalmodule(indexmodule).numberofmontage;
        if sum(indeximagestotalmodule_int) == 0
            paratext6 = 'There are no images to slow because the 2D Jaccard was higher than 0.50 in all the slices.';
            
            para6 = Text(paratext6);
            para6.Style = sectTitle.Style;
            
            add(ch13, para6);
        else            
            firstpageflag = 0; % flag to create the first page of the report with 4 images just one time
            index4images = 1; %index for the images
            
            % Loop for each structure
            for index1 = 1 : size(struct_name_case, 2)
                
                if indeximagestotalmodule(indexmodule).numberofmontage(index1) == 0 % there are no images for this structure
                    
                else                    
                    % Loading the collages
                    %     imagefiles = dir([struct_name_case{index1} '_montage_*.tiff']);
                    %     dirOutput = dir(fullfile(fileFolder,'GTV*.jpg'));
                    
                    % Folder with the images
                    if ispc
                        fileFolder = [handles.FIELDRT.env.userDir '\FIELDRT_Reports\Images\' module '\' struct_name_case{index1} '\'];
                        dirOutput = dir(fullfile(fileFolder, '*_montage_*.tiff'));
                    else
                        fileFolder = [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/'];
                        dirOutput = dir(fullfile(fileFolder, '*_montage_*.tiff'));
                    end                    
                    
                    % Loop for the number of collages
                    for index2 = 1 : indeximagestotalmodule(indexmodule).numberofmontage(index1)
                        
                        if firstpageflag == 0 % just for the first structure
                            
                            if index2 == 1 % just for the first structure to load the 2x2 image
                                
                                paratext7 = 'In this section, slices where the 2D Jaccard was below a threshold (0.50) are showed.';
                                paratext81 = 'Under contoured regions are represented by cyan contours.';
                                paratext82 = 'Over contoured regions are represented by red contours.';
                                paratext83 = 'Regions of agreement are represented by green contours.';
                                
                                para7 = Text(paratext7);
                                para7.Style = sectTitle.Style;
                                
                                para81 = Text(paratext81);
                                para81.Style = sectTitle.Style;
                                
                                para82 = Text(paratext82);
                                para82.Style = sectTitle.Style;
                                
                                para83 = Text(paratext83);
                                para83.Style = sectTitle.Style;
                                
                                % Adding the text to the chapter
                                add(ch13,para7);
                                add(ch13,para81);
                                add(ch13,para82);
                                add(ch13,para83);
                                
                                % title
                                structtitleover = clone(parTitle_image);
                                append(structtitleover, int2str(modulenum));
                                append(structtitleover,'.2.');
                                append(structtitleover, int2str(index4images));
                                append(structtitleover,'. ');
                                append(structtitleover, struct_name_case{index1});
                                sec13 = Section('Title', structtitleover);
                                
                                % Adding a line break
                                add(sec13,l);
                                
                                img = Image([fileFolder dirOutput(index2).name]);
                                img.Style = {ScaleToFit(false), HAlign('center'), Height('5.59in')};
                                % img.Style = {HAlign('center')}; %, Height('13in')};
                                
                                structnumb = Text([struct_name_case{index1} ' - (' int2str(index2) '/' int2str(indeximagestotalmodule(indexmodule).numberofmontage(index1)) ')']);
                                structnumb.Style = sectTitle.Style;
                                structnumb.Style = {VAlign('bottom'), HAlign('right')};
                                
                                % Adding everything
                                add(sec13, img);
                                
                                % Adding line breaks
                                add(ch13,l);
                                add(ch13,l);    
                                
                                add(sec13,structnumb);
                                
                                add(ch13, sec13);
                                
                                index4images = index4images + 1;
                                
                            else
                                
                                sec14 = Section();
                                
                                
                                img = Image([fileFolder dirOutput(index2).name]);
                                img.Style = {ScaleToFit(false), HAlign('center'), Height('8.05in')};
                                % img.Style = {HAlign('center')}; %, Height('16in')};
                                
                                structnumb = Text([struct_name_case{index1} ' - (' int2str(index2) '/' int2str(indeximagestotalmodule(indexmodule).numberofmontage(index1)) ')']);
                                structnumb.Style = sectTitle.Style;
                                structnumb.Style = {VAlign('bottom'), HAlign('right')};
                                
                                % Adding everything
                                add(sec14, img);
                                
                                add(sec14, structnumb);
                                
                                % Adding line breaks
                                add(sec14,l);
                                
                                add(ch13, sec14);
                            end
                        end
                        
                        if firstpageflag == 1 && index4images > 1 % just for the first structure to load the 2x2 image
                            if index2 == 1 %for the first collage
                                
                                % title
                                structtitleover = clone(parTitle_image);
                                append(structtitleover, int2str(modulenum));
                                append(structtitleover,'.2.');
                                append(structtitleover, int2str(index4images));
                                append(structtitleover,'. ');
                                append(structtitleover, struct_name_case{index1});
                                sec15 = Section('Title', structtitleover);
                                
                                
                                img = Image([fileFolder dirOutput(index2).name]);
                                img.Style = {ScaleToFit(false), HAlign('center'), Height('8.05in')};
                                % img.Style = {HAlign('center')}; %, Height('13in')};
                                
                                structnumb = Text([struct_name_case{index1} ' - (' int2str(index2) '/' int2str(indeximagestotalmodule(indexmodule).numberofmontage(index1)) ')']);
                                structnumb.Style = sectTitle.Style;
                                structnumb.Style = {VAlign('bottom'), HAlign('right')};
                                
                                % Adding everything
                                add(sec15, img);
                                add(sec15,structnumb);
                                
                                add(ch13, sec15);
                                
                                index4images = index4images + 1;
                                
                            else
                                
                                sec16 = Section();
                                
                                
                                img = Image([fileFolder dirOutput(index2).name]);
                                img.Style = {ScaleToFit(false), HAlign('center'), Height('8.05in')};
                                % img.Style = {HAlign('center')}; %, Height('16in')};
                                
                                structnumb = Text([struct_name_case{index1} ' - (' int2str(index2) '/' int2str(indeximagestotalmodule(indexmodule).numberofmontage(index1)) ')']);
                                structnumb.Style = sectTitle.Style;
                                structnumb.Style = {VAlign('bottom'), HAlign('right')};
                                
                                % Adding everything
                                add(sec16, img);
                                
                                add(sec16,structnumb);
                                
                                % Adding a line break
                                add(sec16,l);
                                
                                add(ch13, sec16);
                            end
                        end
                    end
                    
                    firstpageflag = 1;
                    
                end
            end
        end
    end
    
    %% User outline + OARs
    
    if indexmodule == 2        
        % Updating the waiting bar
        wb = waitbar((indexmodule / 8), wb);
        
        module = 'OARs';
        modulenum = 2;
        
        %% First page with Module description and what it is going to be showed
        chtitleover = clone(chTitle);
        append(chtitleover, AutoNumber('sect1'));
        append(chtitleover,'. ');
        append(chtitleover,'User outline + OARs');
        ch2 = Chapter('Title', chtitleover);
        
        %         if handles.Case == 11
        %             paratext9 = ['This section includes only the slices in which too much OAR has been in included in the regions outlined by the user. ' ...
        %                 'Reference contour (yellow), user contour (cyan) are showed with the following OARs (when present): '...
        %                 'vertebra (dark green), aorta (red), right lung (orange), pericardium/great vessels (light green), ' ...
        %                 'liver (pink), stomach (green), azygous vein (blue), left main bronchus (pink) and left lung (white).'];
        %         end
        %
        %         if handles.Case == 31
        %             paratext9 = ['This section includes only the slices in which too much OAR has been in included in the regions outlined by the user. ' ...
        %                 'Reference contour (yellow), user contour (cyan) are showed with the following OARs (when present): '...
        %                 'Bladder (dark green), Bowel (red), Lt Femoral Head (orange), Penile bulb (pink), ' ...
        %                 'Rectum (magenta), Rt Femoral Head (blue).'];
        %         end
        
        paratext9_temp = handles.FIELDRTGSCases(handles.Case(1,1)).PDFOARsdescription(handles.Case(1,2)).Description;
        paratext9 = char(paratext9_temp); % converting the cell to a string
        para9 = Text(paratext9);
        para9.Style = sectTitle.Style;
        
        add(ch2, para9);
        
        % Adding a line break
        add(ch2,l);
        
        
        %% Adding the flagged images   
        %         if indexcase == 11 % Oesophagus Case1
        %             struct_name_case = {'GTV', 'CTVB'};
        %         end
        %
        %         if indexcase == 31 % Prostate Case1
        %
        %         end
        
        struct_name_case = [];
        % Retrieving structures for this module
        struct_name_case = handles.FIELDRTGSCases(handles.Case(1,1)).StructOARs(handles.Case(1,2)).Cases;
        
        indeximagestotalmodule_int = indeximagestotalmodule(indexmodule).numberofmontage;
        if sum(indeximagestotalmodule_int) == 0
            paratext6 = 'There are no images to slow because no OARs have been included in the regions outlined by the user.';
            
            para6 = Text(paratext6);
            para6.Style = sectTitle.Style;
            
            add(ch2, para6);
        else
            
            firstpageflag = 0; % flag to create the first page of the report with 4 images just one time
            index4images = 1; %index for the images
            
            % Loop for each structure
            for index1 = 1 : size(struct_name_case, 2)
                
                if indeximagestotalmodule(indexmodule).numberofmontage(index1) == 0 % there are no images for this structure
                    
                else
                    
                    % Loading the collages
                    %     imagefiles = dir([struct_name_case{index1} '_montage_*.tiff']);
                    %     dirOutput = dir(fullfile(fileFolder,'GTV*.jpg'));
                    
                    % Folder with the images
                    if ispc
                        fileFolder = [handles.FIELDRT.env.userDir '\FIELDRT_Reports\Images\' module '\' struct_name_case{index1} '\'];
                        dirOutput = dir(fullfile(fileFolder, '*_montage_*.tiff'));
                    else
                        fileFolder = [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/'];
                        dirOutput = dir(fullfile(fileFolder, '*_montage_*.tiff'));
                    end
                    
                    
                    % Loop for the number of collages
                    for index2 = 1 : indeximagestotalmodule(indexmodule).numberofmontage(index1)
                        
                        if firstpageflag == 0 % just for the first structure
                            
                            if index2 == 1 % just for the first structure to load the 2x2 image
                                % title
                                structtitleover = clone(parTitle);
                                append(structtitleover, int2str(modulenum));
                                append(structtitleover,'.');
                                append(structtitleover, int2str(index4images));
                                append(structtitleover,'. ');
                                append(structtitleover, struct_name_case{index1});
                                sec13 = Section('Title', structtitleover);
                                
                                
                                img = Image([fileFolder dirOutput(index2).name]);
                                img.Style = {ScaleToFit(false), HAlign('center'), Height('5.59in')};
                                % img.Style = {HAlign('center')}; %, Height('13in')};
                                
                                structnumb = Text([struct_name_case{index1} ' - (' int2str(index2) '/' int2str(indeximagestotalmodule(indexmodule).numberofmontage(index1)) ')']);
                                structnumb.Style = sectTitle.Style;
                                structnumb.Style = {VAlign('bottom'), HAlign('right')};
                                
                                % Adding everything
                                add(sec13, img);
                                
                                
                                add(sec13,structnumb);
                                
                                add(ch2, sec13);
                                
                                index4images = index4images + 1;
                                
                            else
                                
                                sec14 = Section();
                                
                                
                                img = Image([fileFolder dirOutput(index2).name]);
                                img.Style = {ScaleToFit(false), HAlign('center'), Height('8.05in')};
                                % img.Style = {HAlign('center')}; %, Height('16in')};
                                
                                structnumb = Text([struct_name_case{index1} ' - (' int2str(index2) '/' int2str(indeximagestotalmodule(indexmodule).numberofmontage(index1)) ')']);
                                structnumb.Style = sectTitle.Style;
                                structnumb.Style = {VAlign('bottom'), HAlign('right')};
                                
                                % Adding everything
                                add(sec14, img);
                                
                                add(sec14, structnumb);
                                
                                % Adding line breaks
                                add(sec14,l);
                                
                                add(ch2, sec14);
                            end
                        end
                        
                        if firstpageflag == 1 && index4images > 1 % just for the first structure to load the 2x2 image
                            if index2 == 1 %for the first collage
                                
                                % title
                                structtitleover = clone(parTitle);
                                append(structtitleover, int2str(modulenum));
                                append(structtitleover,'.');
                                append(structtitleover, int2str(index4images));
                                append(structtitleover,'. ');
                                append(structtitleover, struct_name_case{index1});
                                sec15 = Section('Title', structtitleover);
                                
                                
                                img = Image([fileFolder dirOutput(index2).name]);
                                img.Style = {ScaleToFit(false), HAlign('center'), Height('7.90in')};
                                % img.Style = {HAlign('center')}; %, Height('13in')};
                                
                                structnumb = Text([struct_name_case{index1} ' - (' int2str(index2) '/' int2str(indeximagestotalmodule(indexmodule).numberofmontage(index1)) ')']);
                                structnumb.Style = sectTitle.Style;
                                structnumb.Style = {VAlign('bottom'), HAlign('right')};
                                
                                % Adding everything
                                add(sec15, img);
                                add(sec15,structnumb);
                                
                                add(ch2, sec15);
                                
                                index4images = index4images + 1;
                                
                            else
                                
                                sec16 = Section();
                                
                                
                                img = Image([fileFolder dirOutput(index2).name]);
                                img.Style = {ScaleToFit(false), HAlign('center'), Height('8.05in')};
                                % img.Style = {HAlign('center')}; %, Height('16in')};
                                
                                structnumb = Text([struct_name_case{index1} ' - (' int2str(index2) '/' int2str(indeximagestotalmodule(indexmodule).numberofmontage(index1)) ')']);
                                structnumb.Style = sectTitle.Style;
                                structnumb.Style = {VAlign('bottom'), HAlign('right')};
                                
                                % Adding everything
                                add(sec16, img);
                                
                                add(sec16,structnumb);
                                
                                % Adding a line break
                                add(sec16,l);
                                
                                add(ch2, sec16);
                            end
                        end
                    end
                    
                    firstpageflag = 1;
                    
                end
            end
        end
    end
   
    %% Min and max    
    if indexmodule == 3
        
        % Updating the waiting bar
        wb = waitbar((indexmodule / 8), wb);
        
        module = 'Minmax';
        modulenum = 3;

        %% First page with Module description and what it is going to be showed
        chtitleover = clone(chTitle);
        append(chtitleover, AutoNumber('sect1'));
        append(chtitleover,'. ');
        append(chtitleover,'Min and max');
        ch3 = Chapter('Title', chtitleover);
        
        paratext10 = 'This section shows the user and the maximum/maximum acceptable outlining area.';
        
        % Style of the text
        para10 = Text(paratext10);
        para10.Style = sectTitle.Style;
        
        paratext11 = 'Reference volume (represented as a yellow contour) is what we consider the ideal volume but there is variability.';
        
        % Style of the text
        para11 = Text(paratext11);
        para11.Style = sectTitle.Style;
         
        paratext12 = 'Maximum volume (represented as a red contour) is what we consider the maximum acceptable outlining area.';
        
        % Style of the text
        para12 = Text(paratext12);
        para12.Style = sectTitle.Style;
        
        paratext13 = 'Minimum volume (represented as a green contour) is what we consider the minimum acceptable outlining area.';
        
        % Style of the text
        para13 = Text(paratext13);
        para13.Style = sectTitle.Style;
        
        paratext14 = 'Please have a look at the following image for a visual explanation of this concept.';
        
        % Style of the text
        para14 = Text(paratext14);
        para14.Style = sectTitle.Style;
        
        % Loading the image of the min and max
        if ispc% Folder with the images
            imgminmax = Image([handles.FIELDRT.env.userDir '\FIELDRT-master\Utilities\Min_max_figure.png']);
        else
            imgminmax = Image([handles.FIELDRT.env.userDir '/FIELDRT-master/Utilities/Min_max_figure.png']);
        end
        
        % Defining the style for the image
        imgminmax.Style = {ScaleToFit(false), HAlign('center'), Height('3.8in')};
        % img.Style = {HAlign('center')}; %, Height('13in')};
        
        paratext15 = 'Ideally your volume will be somewhere in between the maximum and minimum area.';
        
        % Style of the text
        para15 = Text(paratext15);
        para15.Style = sectTitle.Style;
        
        % Adding everything to the chapter
        add(ch3, para10);
        
        % Adding a line break
        add(ch3,l);
        
        add(ch3, para11);
        add(ch3, para12);
        add(ch3, para13);
        
        
        % Adding line breaks
        add(ch3,l);
        add(ch3,l);
        
        add(ch3, para14);
        
        % Adding a line break
        add(ch3,l);
        
        add(ch3, imgminmax);
        
        % Adding a line break
        add(ch3,l);
        
        add(ch3, para15);
        
        % Adding line breaks
        add(ch3,l);
        
        %% Adding the flagged images
        
        %         if indexcase == 11 % Oesophagus Case1
        %                 struct_name_case = {'GTV', 'CTVA', 'CTVB', 'CTVC', 'PTV'};
        %         end
        %
        %         if indexcase == 31 % Prostate Case1
        %                 struct_name_case = {'CTVp', 'CTVpsv'};
        %         end
        
        struct_name_case = [];
        % Retrieving structures for this module
        struct_name_case = handles.FIELDRTGSCases(handles.Case(1,1)).StructMinmax(handles.Case(1,2)).Cases;
        
        indeximagestotalmodule_int = indeximagestotalmodule(indexmodule).numberofmontage;
        if sum(indeximagestotalmodule_int) == 0
            paratext61 = 'There are no images to show. To see your performance for the min and max module please refer to the FIELD';
            para61 = Paragraph(paratext61);
            para61.Style = parTitle_minmaxFIELDRT.Style;
            para61.Style = {OutlineLevel(9)}; %to avoid the paragraph to appear in the table of content
            
            paratext62 = 'RT';
            para62 = Text(paratext62);
            para62.Style = parTitle_minmaxFIELDRT.Style; % To make it equal to the rest
            para62.Style = {VerticalAlign('superscript'),FontSize('9pt')}; % To have it as a superscript
            
            paratext63 = ' viewer.';
            para63 = Text(paratext63);
            para63.Style = parTitle_minmaxFIELDRT.Style;   
            
            append(para61, para62);
            append(para61, para63);
            
            add(ch3, para61);
        else
            
            firstpageflag = 0; % flag to create the first page of the report with 4 images just one time --> not needed for the min and max module
            index4images = 1; %index for the images
            
            % Loop for each structure
            for index1 = 1 : size(struct_name_case, 2)
                
                if indeximagestotalmodule(indexmodule).numberofmontage(index1) == 0 % there are no images for this structure
                    
                else
                    
                    % Loading the collages
                    %     imagefiles = dir([struct_name_case{index1} '_montage_*.tiff']);
                    %     dirOutput = dir(fullfile(fileFolder,'GTV*.jpg'));
                    
                    % Folder with the images
                    if ispc
                        fileFolder = [handles.FIELDRT.env.userDir '\FIELDRT_Reports\Images\' module '\' struct_name_case{index1} '\'];
                        dirOutput = dir(fullfile(fileFolder, '*_montage_*.tiff'));
                    else
                        fileFolder = [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/'];
                        dirOutput = dir(fullfile(fileFolder, '*_montage_*.tiff'));
                    end
                    
                    
                    % Loop for the number of collages
                    for index2 = 1 : indeximagestotalmodule(indexmodule).numberofmontage(index1)
                        
                        if firstpageflag == 0 % just for the first structure
                            
                            if index2 == 1 % just for the first structure to load the 2x2 image
                                % title
                                structtitleover = clone(parTitle);
                                append(structtitleover, int2str(modulenum));
                                append(structtitleover,'.');
                                append(structtitleover, int2str(index4images));
                                append(structtitleover,'. ');
                                append(structtitleover, struct_name_case{index1});
                                sec13 = Section('Title', structtitleover);
                                
                                
                                img = Image([fileFolder dirOutput(index2).name]);
                                % img.Style = {ScaleToFit(false), HAlign('center'), Height('5.59in')};
                                img.Style = {ScaleToFit(false), HAlign('center'), Height('7.90in')};
                                % img.Style = {HAlign('center')}; %, Height('13in')};
                                
                                structnumb = Text([struct_name_case{index1} ' - (' int2str(index2) '/' int2str(indeximagestotalmodule(indexmodule).numberofmontage(index1)) ')']);
                                structnumb.Style = sectTitle.Style;
                                structnumb.Style = {VAlign('bottom'), HAlign('right')};
                                
                                % Adding everything
                                add(sec13, img);
                                add(sec13,structnumb);
                                
                                add(ch3, sec13);
                                
                                index4images = index4images + 1;
                                
                            else
                                
                                sec14 = Section();
                                
                                
                                img = Image([fileFolder dirOutput(index2).name]);
                                img.Style = {ScaleToFit(false), HAlign('center'), Height('8.05in')};
                                % img.Style = {HAlign('center')}; %, Height('16in')};
                                
                                structnumb = Text([struct_name_case{index1} ' - (' int2str(index2) '/' int2str(indeximagestotalmodule(indexmodule).numberofmontage(index1)) ')']);
                                structnumb.Style = sectTitle.Style;
                                structnumb.Style = {VAlign('bottom'), HAlign('right')};
                                
                                % Adding everything
                                add(sec14, img);
                                
                                add(sec14, structnumb);
                                
                                % Adding line breaks
                                add(sec14,l);
                                add(sec14,l);
                                
                                add(ch3, sec14);
                            end
                        end
                        
                        if firstpageflag == 1 && index4images > 1 % just for the first structure to load the 2x2 image
                            if index2 == 1 %for the first collage
                                
                                % title
                                structtitleover = clone(parTitle);
                                append(structtitleover, int2str(modulenum));
                                append(structtitleover,'.');
                                append(structtitleover, int2str(index4images));
                                append(structtitleover,'. ');
                                append(structtitleover, struct_name_case{index1});
                                sec15 = Section('Title', structtitleover);
                                
                                
                                img = Image([fileFolder dirOutput(index2).name]);
                                img.Style = {ScaleToFit(false), HAlign('center'), Height('7.90in')};
                                % img.Style = {HAlign('center')}; %, Height('13in')};
                                
                                structnumb = Text([struct_name_case{index1} ' - (' int2str(index2) '/' int2str(indeximagestotalmodule(indexmodule).numberofmontage(index1)) ')']);
                                structnumb.Style = sectTitle.Style;
                                structnumb.Style = {VAlign('bottom'), HAlign('right')};
                                
                                % Adding everything
                                add(sec15, img);
                                add(sec15,structnumb);
                                                                
                                add(ch3, sec15);
                                
                                index4images = index4images + 1;
                                
                            else
                                
                                sec16 = Section();
                                
                                
                                img = Image([fileFolder dirOutput(index2).name]);
                                img.Style = {ScaleToFit(false), HAlign('center'), Height('8.05in')};
                                % img.Style = {HAlign('center')}; %, Height('16in')};
                                
                                structnumb = Text([struct_name_case{index1} ' - (' int2str(index2) '/' int2str(indeximagestotalmodule(indexmodule).numberofmontage(index1)) ')']);
                                structnumb.Style = sectTitle.Style;
                                structnumb.Style = {VAlign('bottom'), HAlign('right')};
                                
                                % Adding everything
                                add(sec16, img);
                                
                                add(sec16,structnumb);
                                
                                % Adding a line break
                                add(sec16,l);
                                
                                add(ch3, sec16);
                            end
                        end
                    end
                    
                    firstpageflag = 1;
                    
                end
            end
        end
    end
end

add(rpt, ch1); % Over/under description page
wb = waitbar(((indexmodule + 1) / 8), wb); % Updating the waiting bar
add(rpt, ch12); % Over/under quantitative evaluation page
add(rpt, ch13); % Over/under qualitative evaluation page
wb = waitbar(((indexmodule + 2) / 8), wb); % Updating the waiting bar
add(rpt, ch2); % OARs module page
wb = waitbar(((indexmodule + 3) / 8), wb); % Updating the waiting bar
add(rpt, ch3); % Min/max module page

wb = waitbar(((indexmodule + 4) / 8), wb);% Updating the waiting bar

% Close the report (required)
close(rpt);

wb = waitbar(((indexmodule + 5) / 8), wb);% Updating the waiting bar

% Close the waiting bar
close(wb)

% Display the report (optional)
if isdeployed == 0 % checking whether the MATLAB or deployed version is running
    % 0 --> MATLAB version
    % 1 --> deployed version
    
    % Display the report (optional)
    rptview(rpt);
else
    if ispc
        winopen([handles.FIELDRT.env.userDir  '\FIELDRT_Reports\' modifiedStr '.pdf']);
    else
        FIELDRTMacopen([handles.FIELDRT.env.userDir  '/FIELDRT_Reports/' modifiedStr '.pdf']);
    end
end

% Deleting the folder with the images
if ispc
    rmdir([handles.FIELDRT.env.userDir '\FIELDRT_Reports\Images'], 's');
else
    rmdir([handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images'], 's');
end

end