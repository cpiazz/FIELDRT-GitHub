function varargout = FIELDRT(varargin)
% FIELDRT MATLAB code for FIELDRT.fig
%      FIELDRT, by itself, creates a new FIELDRT or raises the existing
%      singleton*.
%
%      H = FIELDRT returns the handle to a new FIELDRT or the handle to
%      the existing singleton*.
%
%      FIELDRT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FIELDRT.M with the given input arguments.
%
%      FIELDRT('Property','Value',...) creates a new FIELDRT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before FIELDRT_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to FIELDRT_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above statictext to modify the response to help FIELDRT

% Last Modified by GUIDE v2.5 17-Feb-2021 13:35:37

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
% Concetta Piazzese November 2017

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @FIELDRT_OpeningFcn, ...
                   'gui_OutputFcn',  @FIELDRT_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before FIELDRT is made visible.
function FIELDRT_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to FIELDRT (see VARARGIN)

clc 

% To deactivate the CERRHotKeys and CERRHotKeysRelease callbacks
set(hObject,'KeyReleaseFcn', ' ');
set(hObject,'WindowKeyPressFcn', ' ')

set(hObject,'CloseRequestFcn', {@quitbutton_Callback, handles}) % To avoid closereq Callback  
% Choose default command line output for FIELDRT
handles.output = hObject;

global status path_GS FIELDRT_h

handles.FIELDRT_h = FIELDRT_h;
  
if status == 0 % A warning message appear and the software stops
    h = warndlg('No GS found. FIELDRT software will stop');
    waitfor(h)
    handles.closeFigure = true;
else    
    % handles.listbox.String  = GS_data;
    % handles.listbox.String{1, 1} = '';
    for k = 1 : size(path_GS, 1)
        idx = strfind(path_GS(k).name, '.fieldrt');
        handles.listbox.String{k, 1} = path_GS(k).name(1, 1 : idx - 1);
    end    
end

handles.statictext.String = 'Select one case from the list box';

% To avoid the defeaut selection of the first element of the list
set(handles.listbox, 'Min', 0, 'Max', 2, 'Value', []);

% Loading FIELDRT logo
if ispc
    FIELDRT_logo_path = [handles.FIELDRT_h.env.userDir '\FIELDRT-master\Utilities\'];   
elseif isunix
    FIELDRT_logo_path = [handles.FIELDRT_h.env.userDir '/FIELDRT-master/Utilities/'];
end

% axes(axisHandle)
image(imread([FIELDRT_logo_path 'FIELDRT_logo.png']))
axis off
axis image


% % Update handles structure
guidata(hObject, handles);

% UIWAIT makes FIELDRT wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = FIELDRT_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Clearing Matlab path


% Get default command line output from handles structure
varargout{1} = handles.output;

% 
% 
% if (isfield(handles,'closeFigure') && handles.closeFigure)
%     figure1_CloseRequestFcn(hObject, eventdata, handles)
% end


% - Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: delete(hObject) closes the figure

delete(hObject);

% --- Executes on button press in loadbutton.
function loadbutton_Callback(hObject, eventdata, handles)
% hObject    handle to loadbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% CERRImportDCM4CHE;

global planC DateVector FIELDRT_path_FIELDRT AttemptInformation FIELDRTDataAttempts FIELDRT_path status FIELDRTGSCases

% Load the FIELDRTData.fieldrt
FIELDRT_path_FIELDRT = [FIELDRT_path 'FIELDRTData.fieldrt']; 
load(FIELDRT_path_FIELDRT, '-mat');

size_original = size(planC{1, 4}, 2);
% study containing the structure(s) to insert.
[planC, DateVector, AttemptInformation, status, Datatopass] = FIELDRTinsertStructsAttempt(planC, FIELDRT_path_FIELDRT, FIELDRTDataAttempts, handles, status);
% status = 1 --> first attempt --> the user can proceed with the analysis
% status = 2/3 --> attempt already analyzed --> pdf report and the viewer
% status = 4 --> the user is closing the window --> nothing happens
% are opened
if status == 1
    set(handles.figure1,'CloseRequestFcn', {@quitbutton_Callback, handles}) % To avoid closereq callback
    
    if size_original == size(planC{1, 4}, 2) && status == 2
        handles.loadbutton.Enable = 'off';
        handles.analyseandresultsbutton.Enable = 'off';
        handles.aboutbutton.Enable = 'on';
        handles.helpbutton.Enable = 'on';
        handles.quitbutton.Enable = 'on';
        set(handles.listbox, 'Value', []);
        handles.statictext.String = sprintf('No file selected! Please try again');
    else
        handles.statictext.String = sprintf('Structures successfully loaded');
        set(handles.listbox, 'Min', 0, 'Max', 2, 'Value', []);
        handles.listbox.Enable = 'off';
        handles.loadbutton.Enable = 'off';
        handles.analyseandresultsbutton.Enable = 'on';
        handles.aboutbutton.Enable = 'on';
        handles.helpbutton.Enable = 'on';
        handles.quitbutton.Enable = 'on';
    end
else
    if status == 2 || status == 3
        
        set(handles.figure1,'CloseRequestFcn', 'closereq;')% To avoid quitbutton_Callback
        
        h = msgbox('The report with the results will open after closing the viewer');
        waitfor(h);
        
        waitfor(FIELDRTResultvisualization(Datatopass));
        
        % Opening the pdf report
        % Name of the report to open
        modifiedStr = strrep(FIELDRTDataAttempts(AttemptInformation.index).AttemptName(1, 1 : (end - 4)), ' ',  '_');
        
        if ispc
            winopen([handles.FIELDRT_h.env.userDir  '\FIELDRT_Reports\' modifiedStr '.pdf']);
        else
            FIELDRTMacopen([handles.FIELDRT_h.env.userDir  '/FIELDRT_Reports/' modifiedStr '.pdf']);
        end
        
        % Settings before the first interface is open again
        cd(handles.FIELDRT_h.env.userDir) % changes the current folder to newFolder.
        
        FIELDRT
        
%         if ispc
%                         
%             set(handles.figure1,'CloseRequestFcn', 'closereq;')% To avoid quitbutton_Callback
%             
%             h = msgbox('The report with the results will open after closing the viewer');
%             waitfor(h);
%             
%             waitfor(FIELDRTResultvisualization(Datatopass));
%             
%             % Opening the pdf report
%             
%             % Name of the report to open
%             modifiedStr = strrep(handles.AttemptInformation.FileName(1, 1 : (end - 4)), ' ',  '_');
% 
%             rpt = Report([handles.FIELDRT.env.userDir  '\FIELDRT_Reports\' modifiedStr],'pdf');
% 
%             rptview(rpt);
%             
%             cd(handles.FIELDRT_h.env.userDir) % changes the current folder to newFolder.
%             FIELDRT
%             
%         elseif isunix
%            %
%             set(handles.figure1,'CloseRequestFcn', 'closereq;')% To avoid quitbutton_Callback            
%             
%             h = msgbox('The report with the results will open after closing the viewer');
%             waitfor(h);
%             
%             waitfor(FIELDRTResultvisualization(Datatopass));
%             
%              % Opening the pdf report
%             
%             % Name of the report to open
%             modifiedStr = strrep(handles.AttemptInformation.FileName(1, 1 : (end - 4)), ' ',  '_');
% 
%             rpt = Report([handles.FIELDRT.env.userDir  '/FIELDRT_Reports/' modifiedStr],'pdf');
% 
%             rptview(rpt);
%             
%             
%             % Settings before the first interface is open again
%             cd(handles.FIELDRT_h.env.userDir) % changes the current folder to newFolder.
%             FIELDRT
%       end
    end
    % To uncomment before releasing the software
    if status == 4 
        
    end
end


% Check if the structures have already been analyzed
% planC = insertDCMStructFIELDRT(planC);


% --- Executes on button press in analyseandresultsbutton.
function analyseandresultsbutton_Callback(hObject, eventdata, handles)
% hObject    handle to analyseandresultsbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global Case planC DateVector FIELDRTDataAttempts AttemptInformation FIELDRT_path_FIELDRT FIELDRTGSCases

handles.analyseandresultsbutton.Enable = 'off';

handles.statictext.String = sprintf('Computing the analysis...');  % message to user

% Inserting information in the file that keeps tracks of the loadings

sizeloadedattempts = size(FIELDRTDataAttempts, 2);
if sizeloadedattempts == 1 && isempty(FIELDRTDataAttempts(1).AttemptName)
    FIELDRTDataAttempts(1).AttemptName = AttemptInformation.FileName;
    FIELDRTDataAttempts(1).AttemptPath = AttemptInformation.FiledirName;    
    FIELDRTDataAttempts(1).SOPInstanceUID = AttemptInformation.SOPInstanceUID;
    FIELDRTDataAttempts(1).Date.Year = DateVector(1, 1);
    FIELDRTDataAttempts(1).Date.Month = DateVector(1, 2);
    FIELDRTDataAttempts(1).Date.Day = DateVector(1, 3);
    FIELDRTDataAttempts(1).Time.Hour = DateVector(1, 4);
    FIELDRTDataAttempts(1).Time.Minutes = DateVector(1, 5);
else
    FIELDRTDataAttempts(sizeloadedattempts + 1).AttemptName = AttemptInformation.FileName;
    FIELDRTDataAttempts(sizeloadedattempts + 1).AttemptPath = AttemptInformation.FiledirName;
    FIELDRTDataAttempts(sizeloadedattempts + 1).SOPInstanceUID = AttemptInformation.SOPInstanceUID;
    FIELDRTDataAttempts(sizeloadedattempts + 1).Date.Year = DateVector(1, 1);
    FIELDRTDataAttempts(sizeloadedattempts + 1).Date.Month = DateVector(1, 2);
    FIELDRTDataAttempts(sizeloadedattempts + 1).Date.Day = DateVector(1, 3);
    FIELDRTDataAttempts(sizeloadedattempts + 1).Time.Hour = DateVector(1, 4);
    FIELDRTDataAttempts(sizeloadedattempts + 1).Time.Minutes = DateVector(1, 5);
end

% Num of structures
indexS = planC{end};
structfieldnum = indexS.structures;
numofStruct = size(planC{1, structfieldnum}, 2);

% Analysing Data
[planC, jaccard3Doutput, jaccard2Doutput, volumestat_GS, volumestat_US, COM_GS, COM_US, volume_Ratio, OARs_finalmask, OARsdummy_finalmask] = FIELDRTAnalysis(planC, Case, FIELDRTGSCases);

% % Coverting all the structs in the axial view to sagittal and coronal
% [sagmasks_structs, cormasks_structs] = FIELDRTconvertstructsview(planC);

handles.statictext.String = sprintf('Opening the viewer'); % message to user

pause(1);

% Setting the folder where to save planC
if ispc
    if ~exist([handles.FIELDRT_h.env.userDir '\FIELDRT-master\MatAttempts'] ,'dir')
        mkdir([handles.FIELDRT_h.env.userDir '\FIELDRT-master\MatAttempts']);
    end
else
    if ~exist([handles.FIELDRT_h.env.userDir '/FIELDRT-master/MatAttempts'] ,'dir')
        mkdir([handles.FIELDRT_h.env.userDir '/FIELDRT-master/MatAttempts']);
    end    
end


% Preparing the colours matrix for AORs
% Case = 1; % To be changed based on the case analysed

% % To be changed based on the case analysed
% if Case == 11 % Oesophagus Case1
%     struct_name_case = {'Vertebrae_GS';'Aorta_GS';'Right lung_GS';'Pericardium/great vessels_GS';'Liver_GS';'Stomach_GS';'Azygous vein_GS';'Left main bronchus_GS'; 'Left lung_GS'};
%     colorlab = [0.1  0.5  0.1;            1 0 0;    1 0.6 0.2;             0 1 0;                1 0 1;     0.3 0.8 0;        0 0 1;           1 0.7 0.8;             1 1 1]; % for contour's colour
%                 % forest green,            red,      orange,             green,                magenta,     avocado,         blue,              pink,               white, 
% 
%                 % other colours: gold: [0.9  0.75 0]
%                 % for more colours https://www.rapidtables.com/web/color/RGB_Color.html
%     if size(colorlab, 1) < size(struct_name_case, 1)
%         for indexdum = (size(struct_name_case, 1) -  size(colorlab, 1)) : size(struct_name_case, 1)
%             % Finding a new colour to plot the contour (a
%             % colour different from red, green or blue or from
%             % the radom colours already used to plot other OARs
%             tag = 1;
%             
%             while tag == 1
%                 colosel = rand(1,3);
%                 if (isequal(colosel, [0 0 1]) || isequal(colosel, [0 1 0]) || isequal(colosel, [1 0 0]) || (colosel(1, 1) < 0.2 && colosel(1, 2) < 0.2) ...
%                         || (colosel(1, 2) < 0.2 && colosel(1, 3) < 0.2) || (colosel(1, 1) < 0.2 && colosel(1, 3) < 0.2))
%                     
%                 else
%                     sum = 0;
%                     for indexcolor = 1 : size(colorlab, 1)
%                         if  colosel == colorlab(indexcolor, :)
%                             sum = sum + 1;
%                         else
%                             sum = sum + 0;
%                         end
%                     end
%                     if  sum == 0
%                         colorlab = [colorlab; colosel];
%                         tag = 0;
%                         break
%                     else
%                         
%                     end
%                 end
%                 
%             end
%             
%         end
%     end
% end
% 
% if Case == 31 % Prostate Case1
%     struct_name_case = {'Bladder'; 'Bowel'; 'Lt FemHead'; 'Penile bulb'; 'Rectum'; 'Rt FemHead'};
%     colorlab = [0.1  0.5  0.1;      1 0 0;    1 0.6 0.2;   1 0.7 0.8;     1 0 1;      0 0 1]; % for contour's colour
%                 % forest green,      red,      orange,     pink,         magenta,     blue,            
% 
%                 % other colours: gold: [0.9  0.75 0]
%                 % for more colours https://www.rapidtables.com/web/color/RGB_Color.html
%     if size(colorlab, 1) < size(struct_name_case, 1)
%         for indexdum = (size(struct_name_case, 1) -  size(colorlab, 1)) : size(struct_name_case, 1)
%             % Finding a new colour to plot the contour (a
%             % colour different from red, green or blue or from
%             % the radom colours already used to plot other OARs
%             tag = 1;
%             
%             while tag == 1
%                 colosel = rand(1,3);
%                 if (isequal(colosel, [0 0 1]) || isequal(colosel, [0 1 0]) || isequal(colosel, [1 0 0]) || (colosel(1, 1) < 0.2 && colosel(1, 2) < 0.2) ...
%                         || (colosel(1, 2) < 0.2 && colosel(1, 3) < 0.2) || (colosel(1, 1) < 0.2 && colosel(1, 3) < 0.2))
%                     
%                 else
%                     sum = 0;
%                     for indexcolor = 1 : size(colorlab, 1)
%                         if  colosel == colorlab(indexcolor, :)
%                             sum = sum + 1;
%                         else
%                             sum = sum + 0;
%                         end
%                     end
%                     if  sum == 0
%                         colorlab = [colorlab; colosel];
%                         tag = 0;
%                         break
%                     else
%                         
%                     end
%                 end
%                 
%             end
%             
%         end
%     end
% end

% Preparing all the data required for the visualization panel and the pdf report creation
Datatopass.FIELDRT = handles.FIELDRT_h;
Datatopass.planC = planC;
% Datatopass.sagmasks_structs = sagmasks_structs;
% Datatopass.cormasks_structs = cormasks_structs;
Datatopass.numofStruct = numofStruct; 
Datatopass.jaccard3Doutput = jaccard3Doutput; 
Datatopass.jaccard2Doutput = jaccard2Doutput;
Datatopass.volumestat_GS = volumestat_GS; 
Datatopass.volumestat_US = volumestat_US;
Datatopass.COM_GS = COM_GS; 
Datatopass.COM_US = COM_US; 
Datatopass.volume_Ratio = volume_Ratio; 
% Datatopass.colorlab = colorlab;
Datatopass.OARs_finalmask = OARs_finalmask;
Datatopass.OARsdummy_finalmask = OARsdummy_finalmask;
Datatopass.Case = Case;
Datatopass.FIELDRTDataAttempts = FIELDRTDataAttempts;
Datatopass.AttemptInformation = AttemptInformation;
Datatopass.FIELDRTGSCases = FIELDRTGSCases;

if ispc
    save([handles.FIELDRT_h.env.userDir '\FIELDRT-master\MatAttempts\' FIELDRTDataAttempts(size(FIELDRTDataAttempts, 2)).AttemptName(1 : end - 4) '_' ...
    int2str(FIELDRTDataAttempts(size(FIELDRTDataAttempts, 2)).Date.Day) '_' int2str(FIELDRTDataAttempts(size(FIELDRTDataAttempts, 2)).Date.Month) '_' ...
    int2str(FIELDRTDataAttempts(size(FIELDRTDataAttempts, 2)).Date.Year) '_' int2str(FIELDRTDataAttempts(size(FIELDRTDataAttempts, 2)).Time.Hour) '_' ...
    int2str(FIELDRTDataAttempts(size(FIELDRTDataAttempts, 2)).Time.Minutes) '.mat'], 'Datatopass'); % , '-v7.3'); 
else
    save([handles.FIELDRT_h.env.userDir '/FIELDRT-master/MatAttempts/' FIELDRTDataAttempts(size(FIELDRTDataAttempts, 2)).AttemptName(1 : end - 4) '_' ...
    int2str(FIELDRTDataAttempts(size(FIELDRTDataAttempts, 2)).Date.Day) '_' int2str(FIELDRTDataAttempts(size(FIELDRTDataAttempts, 2)).Date.Month) '_' ...
    int2str(FIELDRTDataAttempts(size(FIELDRTDataAttempts, 2)).Date.Year) '_' int2str(FIELDRTDataAttempts(size(FIELDRTDataAttempts, 2)).Time.Hour) '_' ...
    int2str(FIELDRTDataAttempts(size(FIELDRTDataAttempts, 2)).Time.Minutes) '.mat'], 'Datatopass'); % , '-v7.3'); 
end

% Opening the visualization panel with the resulting contours
if ispc
    cd([handles.FIELDRT_h.env.userDir '\FIELDRT-master\Viewer\']) % changes the current folder to newFolder.
else
    cd([handles.FIELDRT_h.env.userDir '/FIELDRT-master/Viewer/']) % changes the current folder to newFolder. 
end

save(FIELDRT_path_FIELDRT, 'FIELDRTDataAttempts', 'FIELDRTGSCases');

% set(handles.figure1,'CloseRequestFcn', @closereq) % To avoid quitbutton_Callback
% set(handles.figure1,'DeleteFcn', closereq)
% set(handles.figure1,'CloseRequestFcn', '')% To avoid quitbutton_Callback
set(handles.figure1,'CloseRequestFcn', 'closereq;')% To avoid quitbutton_Callback

waitfor(FIELDRTResultvisualization(Datatopass)); 
% 
% cd(handles.FIELDRT_h.env.userDir) % changes the current folder to newFolder.
% 
% Creation of the reduced pdf report (only flagged images for over/under
% and OARs modules, no images for min and max)
FIELDRTPdfreportcreation(Datatopass);

% % Creation of the full pdf report (all images for all modules)
% FIELDRTfullPdfreportcreation(Datatopass);

% % Message that the pdf report has been created in the folder
% 
% 
% % Reopening FIELDRT
FIELDRT
% Enable again the analyseandresults button
% handles.analyseandresultsbutton.Enable = 'on';
% Updating the file that keeps tracks of the loadings
save(FIELDRT_path_FIELDRT, 'FIELDRTDataAttempts', 'FIELDRTGSCases');




% --- Executes on button press in aboutbutton.
function aboutbutton_Callback(hObject, eventdata, handles)
% hObject    handle to aboutbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

hbox = msgbox('Coming soon');
waitfor(hbox);
% web http://cerr.info/cerrwiki/index.php/CERR?w=CERRWiKi


% --- Executes on button press in quitbutton.
function quitbutton_Callback(hObject, eventdata, handles)
% hObject    handle to quitbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% global handles
global FIELDRT_h
% hf = gcf;
% gco(hf); % To retrieve the handles variable

choice = questdlg('Do you want really quit? Unsaved data will be lost', ...
	'Quit', ...
	'Yes','No','No');
% Handle response
switch choice
    case 'Yes'
        
        % clear FIELDRT.m
        closereq;
        
        
        rmpath(genpath(FIELDRT_h.env.userDir))
        
        clc
        clear
        
    case 'No'
end


% --- Executes on button press in helpbutton.
function helpbutton_Callback(hObject, eventdata, handles)
% hObject    handle to helpbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

hbox = msgbox('Coming soon');
waitfor(hbox);


% --- Executes on button press in savebutton.
function savebutton_Callback(hObject, eventdata, handles)
% hObject    handle to savebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB


% --- Executes on button press in viewresultsbutton.
function viewresultsbutton_Callback(hObject, eventdata, handles)
% hObject    handle to viewresultsbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in saveresultsbutton.
function saveresultsbutton_Callback(hObject, eventdata, handles)
% hObject    handle to saveresultsbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in listbox.
function listbox_Callback(hObject, eventdata, handles)
% hObject    handle to listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global path_GS planC FIELDRT_path_FIELDRT Case

% The loadbutton is turned off 
handles.loadbutton.Enable = 'off';

% The user selects a case
index_selected = get(handles.listbox,'Value');
set(handles.listbox, 'Min', 0, 'Max', 1); % to avoid multiple selections

% Deleting planC (if the user select another case) (?????)
planC = [];

% Loading the GS case based on which list box item the user selected

if ispc
    load([path_GS(index_selected).folder '\' path_GS(index_selected).name], '-mat');
elseif isunix
    load([path_GS(index_selected).folder '/' path_GS(index_selected).name], '-mat');
end

% Oesophagus Case1 --> [1 1], Oesophagus Case1 --> [1 2], Oesophagus Case3 --> [1 3], Oesophagus Case4 --> [1 4].
% Pelvis Case1 --> [2 1], Pelvis Case2 --> [2 2], Pelvis Case3 --> [2 3].
% Prostate Case1 --> [3 1], Prostate Case2 --> [3 2], Prostate Case3 --> [3 3]
% if strcmp(handles.listbox.String{index_selected, 1}, 'Oesophagus Case1')
%     Case = [1 1];
% end
% 
% if strcmp(handles.listbox.String{index_selected, 1}, 'Oesophagus Case4')
%     Case = [1 4];
% end

if strcmp(handles.listbox.String{index_selected, 1}, 'Prostate Case1')
    Case = [1 1];
end

% Renaming the structures present in planC
planC = FIELDRTrenameGSstructures(planC);

% Renaming the structures present in planC
idx = strfind(path_GS(index_selected).name, '.fieldrt');
handles.statictext.String = sprintf([path_GS(index_selected).name(1, 1 : (idx - 1)) ' loaded! \n Click the submit button']);  
handles.loadbutton.Enable = 'on';
% planC{1, 3}.uniformScanInfo.scanFileName  

% Hints: contents = cellstr(get(hObject,'String')) returns listbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox


% --- Executes during object creation, after setting all properties.
function listbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

global FIELDRT_path status path_GS FIELDRT_h

% Obtaining the path of the script
fullpath = mfilename('fullpath'); 
cd(fullpath(1 : (end - 7)));

% Load FIELDRT configuration file
FIELDRT_h = FIELDRT_makeConfig;
% handles.FIELDRT_h = FIELDRT_h;

if ispc
    % Just because char(java.lang.System.getProperty('user.dir'))doesn't
    % give back the current working folder. To delete in the future if the 
    % problem is fixed by java developers.
    FIELDRT_h.env.userDir = cd;
    save(FIELDRT_h.FIELDRTConfigFilename,'FIELDRT_h');
elseif isunix
    FIELDRT_h.env.userDir = cd;
    save(FIELDRT_h.FIELDRTConfigFilename,'FIELDRT_h');
end

% % Find FIELDRT.root
% str = mfilename('fullpath');
status = 1;

% Initializing the path. Added here because listbox_CreateFcn function is
% the first to run
if ispc
%     indV = find(str == '\');
%     ind = max(indV);
    FIELDRT_path = [FIELDRT_h.env.userDir '\FIELDRT-master\'];   
elseif isunix
%     indV = find(str == '/');
%     ind = max(indV);
    FIELDRT_path = [FIELDRT_h.env.userDir '/FIELDRT-master/'];
end

if isdeployed == 0 % checking whether the MATLAB or deployed version is running
                            
    % Adding the folder to path
    addpath(genpath(FIELDRT_h.env.userDir))
    addpath(genpath(FIELDRT_path));
end

% [FIELDRT_path] = FIELDRTgetFIELDRTpath(str);

GS_data = [];

% Retrieving GS data for the list
[path_GS] = FIELDRTgetGSData(FIELDRT_path);

% If no GS cases are found
if  size(path_GS, 1) == 0
    status = 0;
else    
    % Hint: listbox controls usually have a white background on Windows.
    %       See ISPC and COMPUTER.
    if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
        set(hObject,'BackgroundColor','white');
    end
end

% handles    structure with handles and user data (see GUIDATA)
