function varargout = FIELDRTResultvisualization(varargin)
% FIELDRTRESULTVISUALIZATION MATLAB code for FIELDRTResultvisualization.fig
%      FIELDRTRESULTVISUALIZATION, by itself, creates a new FIELDRTRESULTVISUALIZATION or raises the existing
%      singleton*.
%
%      H = FIELDRTRESULTVISUALIZATION returns the handle to a new FIELDRTRESULTVISUALIZATION or the handle to
%      the existing singleton*.
%
%      FIELDRTRESULTVISUALIZATION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FIELDRTRESULTVISUALIZATION.M with the given input arguments.
%
%      FIELDRTRESULTVISUALIZATION('Property','Value',...) creates a new FIELDRTRESULTVISUALIZATION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before FIELDRTResultvisualization_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to FIELDRTResultvisualization_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

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
% Concetta Piazzese December 2017

% Edit the above text to modify the response to help FIELDRTResultvisualization

% Last Modified by GUIDE v2.5 17-Feb-2021 12:46:28

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @FIELDRTResultvisualization_OpeningFcn, ...
                   'gui_OutputFcn',  @FIELDRTResultvisualization_OutputFcn, ...
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


% --- Executes just before FIELDRTResultvisualization is made visible.
function FIELDRTResultvisualization_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to FIELDRTResultvisualization (see VARARGIN)

% set(hObject,'CloseRequestFcn', @myclosereq) % To avoid closereq Callback  

% Setting the Pan string (UNICODE is not working if added in the GUI)
[b, map] = imread('hand.gif');
a = ind2rgb(b,map);
[r,c,d]=size(a); 
x=ceil(r/35); 
y=ceil(c/35); 
g=a(1:x:end,1:y:end,:);
g(g==255)=5.5*255;
set(handles.pantogglebutton,'CData',g);

% set(hObject,'CloseRequestFcn', {@quitbutton_Callback, handles}) % To avoid closereq Callback  

% Choose default command line output for FIELDRTResultvisualization
handles.output = hObject;

% set(ll,'CloseRequestFcn', @closereq) % To avoid quitbutton_Callback
% close FIELDRT
global status

if status == 1
    global planC
    
    planC = varargin{1, 1}.planC;
    
    handles.numofStruct = varargin{1, 1}.numofStruct;
    handles.FIELDRT = varargin{1, 1}.FIELDRT;
    handles.jaccard3Doutput = varargin{1, 1}.jaccard3Doutput;
    handles.jaccard2Doutput = varargin{1, 1}.jaccard2Doutput;
    handles.volumestat_GS = varargin{1, 1}.volumestat_GS;
    handles.volumestat_US = varargin{1, 1}.volumestat_US;
    handles.COM_GS = varargin{1, 1}.COM_GS;
    handles.COM_US = varargin{1, 1}.COM_US;
    handles.volume_Ratio = varargin{1, 1}.volume_Ratio;
    % handles.colorlab = varargin{1, 1}.colorlab;
    handles.OARs_finalmask = varargin{1, 1}.OARs_finalmask;
    handles.OARsdummy_finalmask = varargin{1, 1}.OARsdummy_finalmask;
    handles.Case = varargin{1, 1}.Case;
%     handles.sagmasks_structs = varargin{1, 1}.sagmasks_structs;
%     handles.cormasks_structs = varargin{1, 1}.cormasks_structs;
    handles.FIELDRTGSCases = varargin{1, 1}.FIELDRTGSCases;
end

if status == 2 || status == 3
    planC = varargin{1, 1}.planC;
    
    handles.numofStruct = varargin{1, 1}.numofStruct;
    handles.FIELDRT = varargin{1, 1}.FIELDRT;
    handles.jaccard3Doutput = varargin{1, 1}.jaccard3Doutput;
    handles.jaccard2Doutput = varargin{1, 1}.jaccard2Doutput;
    handles.volumestat_GS = varargin{1, 1}.volumestat_GS;
    handles.volumestat_US = varargin{1, 1}.volumestat_US;
    handles.COM_GS = varargin{1, 1}.COM_GS;
    handles.COM_US = varargin{1, 1}.COM_US;
    handles.volume_Ratio = varargin{1, 1}.volume_Ratio;
    % handles.colorlab = varargin{1, 1}.colorlab;
    handles.OARs_finalmask = varargin{1, 1}.OARs_finalmask; 
    handles.OARsdummy_finalmask = varargin{1, 1}.OARsdummy_finalmask;
    handles.Case = varargin{1, 1}.Case;
    handles.FIELDRTGSCases = varargin{1, 1}.FIELDRTGSCases;
end


% set(handles.figure1,'CloseRequestFcn', @closereq) % To avoid quitbutton_Callback


indexS = planC{end};
structfieldnum = indexS.structures;

% to avoid multiple selections
set(handles.listbox1, 'Min', 0, 'Max', 2, 'Value', []);


% Setting the handles for the zoom button
handles.zoomvalue = 1;


% To be changed based on the case analysed
% Case = 1;
% 
% if Case == 1 
%     struct_name_case = {'GTV_GS', 'CTVA_GS', 'CTVB_GS', 'CTVC_GS', 'PTV_GS'};    
% else
%     
% end
% 
% for k = 1 : size(struct_name_case, 2) % handles.numofStruct / 2
%     struct_name_temp = struct_name_case{k};
%     handles.listbox1.String{k, 1} = struct_name_temp(1 : end - 3);
% end
% % To avoid the defeaut selection of the first element of the list
% set(handles.listbox1, 'Min', 0, 'Max', 2, 'Value', []);

% set(handles.figure1,'CloseRequestFcn', @closereq) % To avoid quitbutton_Callback


%% Create the scrollabel figure
scanfieldnum = indexS.scan;
volsize = size(planC{1, scanfieldnum}.scanArray);

% Retrieving some information
indexS = planC{end};
handles.structfieldnum = indexS.structures;
handles.scanfieldnum = indexS.scan;
handles.volsize = size(planC{1, handles.scanfieldnum}.scanArray);
handles.planC = planC;

handles.currentSlice_ax = round(handles.volsize(1, 3) / 2);
set(handles.numslice, 'String', [int2str(handles.currentSlice_ax) '/' int2str(size(handles.planC{1, handles.scanfieldnum}.scanArray, 3))]);

% % Creating the rgb image for filling the plotted contours!
% red = cat(3, ones(volsize(1 : 2)), zeros(volsize(1 : 2)), zeros(volsize(1 : 2)));
% green = cat(3, zeros(volsize(1 : 2)), ones(volsize(1 : 2)), zeros(volsize(1 : 2)));
% blue = cat(3, zeros(volsize(1 : 2)), zeros(volsize(1 : 2)), ones(volsize(1 : 2)));
% handles.hShow = imshow(red, []);
% handles.hShow = imshow(green, []);
% handles.hShow = imshow(blue, []);
hold(handles.axes3, 'off');

% % Range of grey level intensities
% minSlice = min(min(handles.planC{1, handles.scanfieldnum}.scanArray(:, :, handles.currentSlice_ax)));
% maxSlice = max(max(handles.planC{1, handles.scanfieldnum}.scanArray(:, :, handles.currentSlice_ax)));

hold(handles.axes5, 'off');
hold(handles.axes6, 'off');
hold(handles.axes3, 'on');

% Visualizing the image in the axial view
axes(handles.axes3)

% Preparing the grey value range to use
scanSet = 1; % to change if there are more than one scan
%CTLevel     = stateS.optS.CTLevel + CTOffset;
%CTWidth     = stateS.optS.CTWidth;
%CTLow       = CTLevel - CTWidth/2;
%CTHigh      = CTLevel + CTWidth/2;
% scanUID = ['c',repSpaceHyp(planC{scanfieldnum}(scanSet).scanUID(max(1,end-61):end))];
CTOffset    = handles.planC{indexS.scan}(scanSet).scanInfo(1).CTOffset;
CTLevel_image = handles.planC{handles.scanfieldnum}(scanSet).scanInfo(1).DICOMHeaders.WindowCenter(end);
CTWidth_image = handles.planC{handles.scanfieldnum}(scanSet).scanInfo(1).DICOMHeaders.WindowWidth(end);
CTLevel     = CTLevel_image + CTOffset;
CTWidth     = CTWidth_image;
CTLow       = CTLevel - CTWidth/2;
CTHigh      = CTLevel + CTWidth/2;
scanMin = single(min(planC{handles.scanfieldnum}(scanSet).scanArray(:)));
scanMax = single(max(planC{handles.scanfieldnum}(scanSet).scanArray(:)));
handles.CTLow = max(CTLow,scanMin);
handles.CTHigh = min(CTHigh,scanMax);
handles.clippedCT_trans = clip(handles.planC{1, handles.scanfieldnum}.scanArray, CTLow, CTHigh, 'limits');

% To avoid the movement of the axes
hold(handles.axes3, 'off');

set(handles.axes3,'Units','pixels');
resizePos_trans = get(handles.axes3,'Position');

% Resizing the image
handles.myImage_trans= imresize(handles.clippedCT_trans(:, :, handles.currentSlice_ax), [resizePos_trans(3) resizePos_trans(4)]);

handles.hShow_trans = imshow(handles.myImage_trans, [handles.CTLow handles.CTHigh], 'Parent', handles.axes3);
% Setting some handles - hShow should always be the last handle to be set!!
% handles.numslice = numslice;
set(handles.hShow_trans,'cdata', handles.clippedCT_trans(:, :, handles.currentSlice_ax));
% Setting the handle for the callback when the mouse wheel is scrolled
% set(gcf, 'WindowScrollWheelFcn', {@FIELDRTwheel,handles});
% % set(handles.axes3,'Units','normalized');

axis off
axis image
axis normal

% Visualizing the image in the sagittal view
axes(handles.axes5)
set(handles.axes5,'Units','pixels');
handles.resizePos_sag = get(handles.axes5,'Position');

handles.ImgSg = squeeze(permute(handles.planC{1, handles.scanfieldnum}.scanArray, [3 1 2 4])); % Sagittal view image
% minSlice_ImgSg = min(min(ImgSg(:, :, handles.currentSlice_ax))); % Range of grey level intensities
% maxSlice_ImgSg = max(max(ImgSg(:, :, handles.currentSlice_ax))); % Range of grey level intensities
handles.currentSlice_sag = round(handles.volsize(1, 2) / 2);
set(handles.numslice2, 'String', [int2str(handles.currentSlice_sag) '/' int2str(size(handles.ImgSg, 3))]);

handles.resizePos_sag(3) = (handles.resizePos_sag(3) / size(handles.ImgSg, 2)) * size(handles.ImgSg, 1);

% Resizing the image
handles.myImage_sag= imresize(handles.ImgSg(:, :, handles.currentSlice_sag), [handles.resizePos_sag(3) handles.resizePos_sag(4)]);

% % To avoid the movement of the axes
% hold(handles.axes5, 'off');
% cla reset

handles.hShow_sag = imshow(handles.myImage_sag, [handles.CTLow  handles.CTHigh], 'Parent', handles.axes5);
% daspect([15 5 1])% aspect ratio for the image
% set(handles.axes5,'Units','normalized');


axis off
axis image
axis normal

% Visualizing the image in the coronal view
axes(handles.axes6)
set(handles.axes6,'Units','pixels');
handles.resizePos_cor = get(handles.axes6,'Position');

% position = get(handles.axes6,'Position');

handles.ImgCr = squeeze(permute(handles.planC{1, handles.scanfieldnum}.scanArray, [3 2 1 4])); % Coronal view image
% minSlice_ImgCr = min(min(ImgCr(:, :, handles.currentSlice_ax))); % Range of grey level intensities
% maxSlice_ImgCr = max(max(ImgCr(:, :, handles.currentSlice_ax))); % Range of grey level intensities
handles.currentSlice_cor = round(handles.volsize(1, 1) / 2);
set(handles.numslice3, 'String', [int2str(handles.currentSlice_cor) '/' int2str(size(handles.ImgCr, 3))]);

% To find the ratio scale of the 
handles.resizePos_sag(3) = (handles.resizePos_cor(3) / size(handles.ImgCr, 2)) * size(handles.ImgCr, 1);

% Resizing the image
handles.myImage_cor= imresize(double(handles.ImgCr(:, :, handles.currentSlice_cor)), [handles.resizePos_cor(3) handles.resizePos_cor(4)]);

% % To avoid the movement of the axes
% hold(handles.axes6, 'off');
% cla reset

handles.hShow_cor = imshow(handles.myImage_cor, [handles.CTLow  handles.CTHigh], 'Parent', handles.axes6);

axis off
axis image
axis normal

% daspect([15 5 1])% aspect ratio for the image
% set(handles.axes6,'Units','normalized');
% set(handles.axes6,'Position', position);
% Reactivating the axial view axes
axes(handles.axes3)  

% handles.hShow  = imshow(handles.axes3.Children.CData);

% handles.hShow  = imshow(handles.planC{1, handles.scanfieldnum}.scanArray(:, :, handles.currentSlice_ax), [minSlice maxSlice]);

% handles.hShow = hShow;

% Initializing the plot handle
index_selected = get(handles.listbox1,'Value');

% Setting the color matrix for plotted contours
% col = [0 0 1; 1 0 0; 0 1 0];
% if isempty(index_selected) 
%     % So to have 3 plot handles
%     for indexstruc = 1 : 3
%         xinit = 0.00001;
%         yinit = 0.00001;
%         handles.hPlot = plot(xinit, yinit, 'Color', col(indexstruc, :));
%         hold on
%         %     set(handles.hPlot, 'XData', xinit, 'YData', yinit) % , 'linewidth', 4,...
%         %                                     % 'color', [0 0 1], 'visible','on')
%         %      hold(handles.axes3, 'on');
%     end
% else 
%     FIELDRTplotContours(handles, index_selected);
% end
FIELDRTplotContours(handles, index_selected);

axes(handles.axes3) 
axis off
axis image
axis normal

% set(handles.axes3,'Units','normalized');

% Retrieving axes limits for the zoom function
handles.axialxlim(1, 1) = handles.axes3.XLim(1,1);
handles.axialxlim(1, 2) = handles.axes3.XLim(1,2);
handles.axialylim(1, 1) = handles.axes3.YLim(1,1);
handles.axialylim(1, 2) = handles.axes3.YLim(1,2);
handles.sagittalxlim(1, 1) = handles.axes5.XLim(1,1);
handles.sagittalxlim(1, 2) = handles.axes5.XLim(1,2);
handles.sagittalylim(1, 1) = handles.axes5.YLim(1,1);
handles.sagittalylim(1, 2) = handles.axes5.YLim(1,2);
handles.coronalxlim(1, 1) = handles.axes6.XLim(1,1);
handles.coronalxlim(1, 2) = handles.axes6.XLim(1,2);
handles.coronalylim(1, 1) = handles.axes6.YLim(1,1);
handles.coronalylim(1, 2) = handles.axes6.YLim(1,2);

% The weel callback has to be defined after you have set all the required
% variables
set(gcf, 'WindowScrollWheelFcn', {@FIELDRTwheel,handles});

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes FIELDRTResultvisualization wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = FIELDRTResultvisualization_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Get default command line output from handles structure
varargout{1} = handles.output;
% 
% function myclosereq
% %CLOSEREQ  Figure close request function.
% %   CLOSEREQ deletes the current figure window.  By default, CLOSEREQ is
% %   the CloseRequestFcn for new figures.
% 
% %   Copyright 1984-2012 The MathWorks, Inc.
% 
% %   Note that closereq now honors the user's ShowHiddenHandles setting
% %   during figure deletion.  This means that deletion listeners and
% %   DeleteFcns will now operate in an environment which is not guaranteed
% %   to show hidden handles.
% if isempty(gcbf)
%     if length(dbstack) == 1
%         warning(message('MATLAB:closereq:ObsoleteUsage'));
%     end
%     close('force');
% else
%     delete(gcbf);
% end
% 
% cd(handles.FIELDRT_h.env.userDir) % changes the current folder to newFolder.
% FIELDRT


% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1

% Pan off
set(handles.pantogglebutton, 'Value', 0);
pan off

% The user selects a case
index_selected = get(handles.listbox1,'Value');
set(handles.listbox1, 'Min', 0, 'Max', 1); % to avoid multiple selections

% Getting the middle slice position of the selected structure
[rV,cV,sV] = xyztom(handles.COM_GS{index_selected, 1}.COM(1, 1),handles.COM_GS{index_selected, 1}.COM(1, 2),handles.COM_GS{index_selected, 1}.COM(1, 3), 1, handles.planC);

handles.currentSlice_ax = round(sV); 

% Plotting the right slice after deleting any possible line object
% Deleting all annotations in the image
axesHandlesToChildObjects = findobj(gca, 'Type', 'Line');
if ~isempty(axesHandlesToChildObjects)
    delete(axesHandlesToChildObjects);
end

set(handles.axes3.Children, 'cdata', handles.planC{1, handles.scanfieldnum}.scanArray(:, :, handles.currentSlice_ax));
set(handles.numslice, 'String', [int2str(handles.currentSlice_ax) '/' int2str(size(handles.planC{1, handles.scanfieldnum}.scanArray, 3))]);

refreshdata(handles.axes3)

FIELDRTwheel(hObject, [], handles);

% Panel with stats (to run just for the Over/under module)
if get(handles.overunder, 'Value') == 1 || get(handles.refvolume, 'Value') == 1
    
    % Volumes
    volGSvalue = num2str(handles.volumestat_GS{index_selected, 1}.vol);
    
    if size(volGSvalue, 2) < 6 % if the number is too long!
        
    else
        volGSvalue = volGSvalue(1 : 6);
    end
    
    volUSvalue = num2str(handles.volumestat_US{index_selected, 1}.vol);
    
    if size(volUSvalue, 2) < 6 % if the number is too long!
        
    else
        volUSvalue = volUSvalue(1 : 6);
    end
    
    % We don't wait this information
    % % Volume Ratio
    % volratiovalue = num2str(handles.volume_Ratio{index_selected, 1});
    %
    % if size(volratiovalue, 2) < 6 % if the number is too long!
    %
    % else
    %     volratiovalue = volratiovalue(1 : 6);
    % end
    %
    % % COM GS used to compute the middle slice of the selected structure
    % COMGSvalue = num2str([num2str(handles.volumestat_GS{index_selected, 1}.COM(1, 1)) ' ' ...
    %               num2str(handles.volumestat_GS{index_selected, 1}.COM(1, 2)) ' ' ...
    %               num2str(handles.volumestat_GS{index_selected, 1}.COM(1, 3))]);
    %
    % % COM US
    % COMUSvalue = num2str([num2str(handles.volumestat_US{index_selected, 1}.COM(1, 1)) ' ' ...
    %               num2str(handles.volumestat_US{index_selected, 1}.COM(1, 2)) ' ' ...
    %               num2str(handles.volumestat_US{index_selected, 1}.COM(1, 3))]);
    
    % 3D Jaccard similarity coefficient
    jaccard3Dvalue = num2str(handles.jaccard3Doutput{index_selected, 1});
    
    if size(jaccard3Dvalue, 2) < 6 % if the number is too long!
        
    else
        jaccard3Dvalue = jaccard3Dvalue(1 : 6);
    end
    
    % 2D Jaccard similarity coefficient
    jaccard2Dvalue = num2str(handles.jaccard2Doutput{index_selected, handles.currentSlice_ax});
    
    if size(jaccard2Dvalue, 2) < 6 % if the number is too long!
        
    else
        jaccard2Dvalue = jaccard2Dvalue(1 : 6);
    end
    
    % Storing all this information in one variable
    quantitatfeed = sprintf('%s\n%s\n%s\n%s\n%s\n%s\n%s\n',['Reference volume [ml] = ' volGSvalue], ...
        ['User volume [ml] = ' volUSvalue], ...
        ['3D Jaccard similarity coefficient = ' jaccard3Dvalue], ...
        ['2D Jaccard similarity coefficient = ' jaccard2Dvalue]);
    %        ['Volume Ratio = ' volratiovalue], ['COM GS = ' COMGSvalue], ['COM US = ' COMUSvalue], ...
    
    set(handles.stats, 'String', quantitatfeed);
end

% handles.oars.Enable = 'on';
% handles.accreg.Enable = 'on';
FIELDRTplotContours(handles, index_selected)

set(gcf, 'WindowScrollWheelFcn', {@FIELDRTwheel,handles});
    
guidata(hObject,handles); % To update the value of the axial axis (you need to use this!!)


% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

global planc

% Closing FIELDRT 
% hf_FIELDRT = findobj('Name','FIELDRT');
% close(hf_FIELDRT)
% hf = findobj('Name', 'FIELDRT');
% set(hf,'DeleteFcn', @closereq)
% delete(hf)
close FIELDRT
 
% ll = gcf;
% set(ll,'CloseRequestFcn', @closereq) % To avoid quitbutton_Callback

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function numslice_Callback(hObject, eventdata, handles)
% hObject    handle to numslice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of numslice as text
%        str2double(get(hObject,'String')) returns contents of numslice as a double

% Pan off
set(handles.pantogglebutton, 'Value', 0);
pan off

% --- Executes during object creation, after setting all properties.
function numslice_CreateFcn(hObject, eventdata, handles)
% hObject    handle to numslice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in zoomresetbutton.
function zoomresetbutton_Callback(hObject, eventdata, handles)
% hObject    handle to zoomresetbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Pan off
set(handles.pantogglebutton, 'Value', 0);
pan off 
    
handles.axes3.XLim(1,1) = handles.axialxlim(1, 1);
handles.axes3.XLim(1,2) = handles.axialxlim(1, 2);
handles.axes3.YLim(1,1) = handles.axialylim(1, 1);
handles.axes3.YLim(1,2) = handles.axialylim(1, 2);

set(handles.Slider_zoom, 'Value', 1);

axesHandlesToChildObjects = findobj(gca, 'Type', 'Line');
if ~isempty(axesHandlesToChildObjects)
    for indexline = 1 : size(axesHandlesToChildObjects, 1)
        if strcmp(axesHandlesToChildObjects(indexline).LineStyle, 'none')
            set(axesHandlesToChildObjects, 'MarkerSize', round(get(handles.Slider_zoom, 'Value'))* 10);
        else
            set(axesHandlesToChildObjects, 'linewidth', round(get(handles.Slider_zoom, 'Value')));
        end
        % set(axesHandlesToChildObjects, 'linewidth', round(get(handles.Slider_zoom, 'Value')));
    end
end

% --- Executes on button press in zoomminusbutton.
function zoomminusbutton_Callback(hObject, eventdata, handles)
% hObject    handle to zoomminusbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Pan off
set(handles.pantogglebutton, 'Value', 0);
pan off

zoom(0.5)

% --- Executes on button press in zoomplusbutton.
function zoomplusbutton_Callback(hObject, eventdata, handles)
% hObject    handle to zoomplusbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Pan off
set(handles.pantogglebutton, 'Value', 0);
pan off

zoom(1.5)

% --- Executes on button press in accreg.
function accreg_Callback(hObject, eventdata, handles) % Min and max
% hObject    handle to accreg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Deselecting all the others toogle buttons
set(handles.oars , 'Value', 0);
set(handles.refvolume , 'Value', 0);
set(handles.overunder , 'Value', 0);

% Pan off
set(handles.pantogglebutton, 'Value', 0);
pan off

% Hint: get(hObject,'Value') returns toggle state of accreg
if handles.accreg.Value == 0 % when toogle is deselected
    
    % to avoid multiple selections
    set(handles.listbox1, 'Min', 0, 'Max', 2, 'Value', []);
    
%     for k = 1 : size(struct_name_case, 2) % handles.numofStruct / 2
%         struct_name_temp = struct_name_case{k};
%         handles.listbox1.String{k, 1} = struct_name_temp(1 : end - 3);
%     end    
    set(handles.listbox1, 'Value', 0);
    handles.listbox1.Enable = 'off';
    index_selected = [];
    % index_selected = get(handles.listbox1, 'Value');
    
    % Deactivating the reference volume control
    set(handles.refvolume_checkbox, 'Value', 0);
    set(handles.stats, 'Value', 0);
    
    % Deactivating the stats panel
    set(handles.stats, 'String', [])
    
    FIELDRTplotContours(handles, index_selected)
else
    
    % Enable the list with the structures
    handles.listbox1.Enable = 'on';
%     
%     if handles.Case == 11
%         struct_name_case = {'GTV_GS', 'CTVA_GS', 'CTVB_GS', 'CTVC_GS', 'PTV_GS'};
%     end
%     
%     if handles.Case == 31
%         struct_name_case = {'CTVp_GS', 'CTVpsv_GS'};
%     end
    
    struct_name_case = handles.FIELDRTGSCases(handles.Case(1,1)).StructMinmax(handles.Case(1,2)).Cases;
    
    for k = 1 : size(struct_name_case, 2) % handles.numofStruct / 2
        %         struct_name_temp = struct_name_case{k};
        %         handles.listbox1.String{k, 1} = struct_name_temp(1 : end - 3);
        handles.listbox1.String{k, 1} = struct_name_case{k};
    end
    
    % to avoid multiple selections
    set(handles.listbox1, 'Min', 0, 'Max', 2, 'Value', []);
    
    % set(handles.listbox1, 'Value', []);
    
    
    index_selected = get(handles.listbox1,'Value');
    
    % Deactivating the stats panel
    set(handles.stats, 'String', [])
    
    % Activating the reference volume control
    set(handles.refvolume_checkbox, 'Value', 1);
    
    % Hint: get(hObject,'Value') returns toggle state of oars
    FIELDRTplotContours(handles, index_selected)
end


% --- Executes on button press in oars.
function oars_Callback(hObject, eventdata, handles) % OARs
% hObject    handle to oars (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Deselecting all the others toogle buttons
set(handles.refvolume , 'Value', 0);
set(handles.accreg , 'Value', 0);
set(handles.overunder , 'Value', 0);

% Deactivating the stats panel
set(handles.stats, 'String', [])

% Pan off
set(handles.pantogglebutton, 'Value', 0);
pan off

if handles.oars.Value == 0 % when toogle is deselected
%     % To be changed based on the case analysed
%      Case = 1;
    
%     if Case == 1
%         struct_name_case = {};
%     else
%         
%     end
    
    % to avoid multiple selections
    set(handles.listbox1, 'Min', 0, 'Max', 2, 'Value', []);
    
%     for k = 1 : size(struct_name_case, 2) % handles.numofStruct / 2
%         struct_name_temp = struct_name_case{k};
%         handles.listbox1.String{k, 1} = struct_name_temp(1 : end - 3);
%     end
    set(handles.listbox1, 'Value', 0);
    handles.listbox1.Enable = 'off';
    index_selected = [];
    % index_selected = get(handles.listbox1, 'Value');
    
    % Deactivating the reference volume control
    set(handles.refvolume_checkbox, 'Value', 0);
    
    FIELDRTplotContours(handles, index_selected)
else
    
    % Enable the list with the structures
    handles.listbox1.Enable = 'on';

%     if handles.Case == 11
%         struct_name_case = {'GTV_GS', 'CTVB_GS'};
%     end 
%     
%     if handles.Case == 31
%         struct_name_case = {'CTVp_GS', 'CTVpsv_GS'};
%     end  

    handles.listbox1.String = [];
    
    struct_name_case = handles.FIELDRTGSCases(handles.Case(1,1)).StructOARs(handles.Case(1,2)).Cases;
    
    for k = 1 : size(struct_name_case, 2) % handles.numofStruct / 2
        %         struct_name_temp = struct_name_case{k};
        %         handles.listbox1.String{k, 1} = struct_name_temp(1 : end - 3);
        handles.listbox1.String{k, 1} = struct_name_case{k};
    end
    
    % to avoid multiple selections
    set(handles.listbox1, 'Min', 0, 'Max', 2, 'Value', []);
    % set(handles.listbox1, 'Value', 1);
    
    index_selected = get(handles.listbox1, 'Value');
    
    % Updating the stats panel
    quantitatfeed = sprintf('',[]);
    
    set(handles.stats, 'String', quantitatfeed);
    
    % Activating the reference volume control
    set(handles.refvolume_checkbox, 'Value', 1);
   
    % Hint: get(hObject,'Value') returns toggle state of oars
    FIELDRTplotContours(handles, index_selected)
end


% --- Executes on button press in overunder.
function overunder_Callback(hObject, eventdata, handles) %Over-under
% hObject    handle to overunder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Deselecting all the others toogle buttons
set(handles.oars , 'Value', 0);
set(handles.accreg , 'Value', 0);
set(handles.refvolume , 'Value', 0);

% Deactivating the stats panel
set(handles.stats, 'String', [])

% Pan off
set(handles.pantogglebutton, 'Value', 0);
pan off

% Hint: get(hObject,'Value') returns toggle state of overunder
if handles.overunder.Value == 0 % when toogle is deselected
%     % To be changed based on the case analysed
%      Case = 1;
    
%     if Case == 1
%         struct_name_case = {};
%     else
%         
%     end
    
    % to avoid multiple selections
    set(handles.listbox1, 'Min', 0, 'Max', 2, 'Value', []);
    
%     for k = 1 : size(struct_name_case, 2) % handles.numofStruct / 2
%         struct_name_temp = struct_name_case{k};
%         handles.listbox1.String{k, 1} = struct_name_temp(1 : end - 3);
%     end    
    set(handles.listbox1, 'Value', 0);
    handles.listbox1.Enable = 'off';    
    index_selected = [];
    % index_selected = get(handles.listbox1, 'Value');
    
    % Deactivating the reference volume control
    set(handles.refvolume_checkbox, 'Value', 0);
    
    FIELDRTplotContours(handles, index_selected)
else
    
    % Enable the list with the structures
    handles.listbox1.Enable = 'on';
    
%     if handles.Case == 11
%         struct_name_case = {'GTV_GS', 'CTVA_GS', 'CTVB_GS', 'CTVC_GS', 'PTV_GS'};
%     end 
%     
%     if handles.Case == 31
%         struct_name_case = {'CTVp_GS', 'CTVpsv_GS'};
%     end  
    
    struct_name_case = handles.FIELDRTGSCases(handles.Case(1,1)).StructOverunder(handles.Case(1,2)).Cases;
    
    for k = 1 : size(struct_name_case, 2) % handles.numofStruct / 2
        %         struct_name_temp = struct_name_case{k};
        %         handles.listbox1.String{k, 1} = struct_name_temp(1 : end - 3);
        handles.listbox1.String{k, 1} = struct_name_case{k};
    end
    
    % to avoid multiple selections
    set(handles.listbox1, 'Min', 0, 'Max', 2, 'Value', []);
    % set(handles.listbox1, 'Value', 1);
    
    index_selected = get(handles.listbox1,'Value');
    
    % Preparing stats feedback
    % Volumes
    if isempty(index_selected)
        
    else
        
        volGSvalue = num2str(handles.volumestat_GS{index_selected, 1}.vol);
        
        if size(volGSvalue, 2) < 6 % if the number is too long!
            
        else
            volGSvalue = volGSvalue(1 : 6);
        end
        
        volUSvalue = num2str(handles.volumestat_US{index_selected, 1}.vol);
        
        if size(volUSvalue, 2) < 6 % if the number is too long!
            
        else
            volUSvalue = volUSvalue(1 : 6);
        end
        
        % 3D Jaccard similarity coefficient
        jaccard3Dvalue = num2str(handles.jaccard3Doutput{index_selected, 1});
        
        if size(jaccard3Dvalue, 2) < 6 % if the number is too long!
            
        else
            jaccard3Dvalue = jaccard3Dvalue(1 : 6);
        end
        
        % 2D Jaccard similarity coefficient
        jaccard2Dvalue = num2str(handles.jaccard2Doutput{index_selected, handles.currentSlice_ax});
        
        if size(jaccard2Dvalue, 2) < 6 % if the number is too long!
            
        else
            jaccard2Dvalue = jaccard2Dvalue(1 : 6);
        end
        
        % Storing all this information in one variable
        quantitatfeed = sprintf('%s\n%s\n%s\n%s\n%s\n%s\n%s\n', ['Reference volume [ml] = ' volGSvalue], ...
            ['User volume [ml] = ' volUSvalue], ...
            ['3D Jaccard similarity coefficient = ' jaccard3Dvalue], ...
            ['2D Jaccard similarity coefficient = ' jaccard2Dvalue]);
        %        ['Volume Ratio = ' volratiovalue], ['COM GS = ' COMGSvalue], ['COM US = ' COMUSvalue], ...
        
        set(handles.stats, 'String', quantitatfeed);
    end
    
    % Activating the reference volume control
    set(handles.refvolume_checkbox, 'Value', 1);
    
    % Hint: get(hObject,'Value') returns toggle state of oars
    FIELDRTplotContours(handles, index_selected)
end


% --- Executes on button press in refvolume.
function refvolume_Callback(hObject, eventdata, handles)
% hObject    handle to refvolume (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Deselecting all the others toogle buttons
set(handles.oars , 'Value', 0);
set(handles.accreg , 'Value', 0);
set(handles.overunder , 'Value', 0);

% Pan off
set(handles.pantogglebutton, 'Value', 0);
pan off

% Hint: get(hObject,'Value') returns toggle state of refvolume
if handles.refvolume.Value == 0 % when toogle is deselected
%     % To be changed based on the case analysed
%      Case = 1;
    
%     if Case == 1
%         struct_name_case = {};
%     else
%         
%     end
    
    % to avoid multiple selections
    set(handles.listbox1, 'Min', 0, 'Max', 2, 'Value', []);
    
    %     for k = 1 : size(struct_name_case, 2) % handles.numofStruct / 2
    %         struct_name_temp = struct_name_case{k};
    %         handles.listbox1.String{k, 1} = struct_name_temp(1 : end - 3);
    %     end
    set(handles.listbox1, 'Value', 0);
    handles.listbox1.Enable = 'off';
    index_selected = [];
    % index_selected = get(handles.listbox1, 'Value');
    
    % Deactivating the reference volume control
    set(handles.refvolume_checkbox, 'Value', 0);
    
    % Deactivating the stats panel
    set(handles.stats, 'String', [])
    
    FIELDRTplotContours(handles, index_selected)
else
    
    % Enable the list with the structures
    handles.listbox1.Enable = 'on';
    
    %     if handles.Case == 11
    %         struct_name_case = {'GTV_GS', 'CTVA_GS', 'CTVB_GS', 'CTVC_GS', 'PTV_GS'};
    %     end
    %
    %     if handles.Case == 31
    %         struct_name_case = {'CTVp_GS', 'CTVpsv_GS'};
    %     end
    
    struct_name_case = handles.FIELDRTGSCases(handles.Case(1,1)).StructOverunder(handles.Case(1,2)).Cases;
    
    for k = 1 : size(struct_name_case, 2) % handles.numofStruct / 2
        %         struct_name_temp = struct_name_case{k};
        %         handles.listbox1.String{k, 1} = struct_name_temp(1 : end - 3);
        handles.listbox1.String{k, 1} = struct_name_case{k};
    end
    
    % to avoid multiple selections
    set(handles.listbox1, 'Min', 0, 'Max', 2, 'Value', []);
    % set(handles.listbox1, 'Value', 1);
    
    index_selected = get(handles.listbox1,'Value');
    
    if isempty(index_selected)
        
    else
        % Preparing stats feedback
        % Volumes
        volGSvalue = num2str(handles.volumestat_GS{index_selected, 1}.vol);
        
        if size(volGSvalue, 2) < 6 % if the number is too long!
            
        else
            volGSvalue = volGSvalue(1 : 6);
        end
        
        volUSvalue = num2str(handles.volumestat_US{index_selected, 1}.vol);
        
        if size(volUSvalue, 2) < 6 % if the number is too long!
            
        else
            volUSvalue = volUSvalue(1 : 6);
        end
        
        % 3D Jaccard similarity coefficient
        jaccard3Dvalue = num2str(handles.jaccard3Doutput{index_selected, 1});
        
        if size(jaccard3Dvalue, 2) < 6 % if the number is too long!
            
        else
            jaccard3Dvalue = jaccard3Dvalue(1 : 6);
        end
        
        % 2D Jaccard similarity coefficient
        jaccard2Dvalue = num2str(handles.jaccard2Doutput{index_selected, handles.currentSlice_ax});
        
        if size(jaccard2Dvalue, 2) < 6 % if the number is too long!
            
        else
            jaccard2Dvalue = jaccard2Dvalue(1 : 6);
        end
        
        % Deactivating the stats panel
        set(handles.stats, 'String', [])
        
        % Activating the reference volume control
        set(handles.refvolume_checkbox, 'Value', 1);
    end
    
    % Activating the reference volume control
    set(handles.refvolume_checkbox, 'Value', 1);
    
    % Hint: get(hObject,'Value') returns toggle state of oars
    FIELDRTplotContours(handles, index_selected)
end

% --- Executes on button press in pantogglebutton.
function pantogglebutton_Callback(hObject, eventdata, handles)
% hObject    handle to pantogglebutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of pantogglebutton
axes(handles.axes3)

% Pan on
if get(handles.pantogglebutton, 'Value') == 1
    pan on
end

% Pan off
if get(handles.pantogglebutton, 'Value') == 0
    pan off
end


% --- Executes on slider movement.
function Slider_zoom_Callback(hObject, eventdata, handles)
% hObject    handle to Slider_zoom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
set(handles.pantogglebutton, 'Value', 0);
pan off

zoom_value = get(handles.Slider_zoom, 'Value');

% WE ARE NOT GOING TO USE THE ZOOM FUNCTION BECAUSE IT APPLIES A RELATIVE
% ZOOMING AND NOT AND ABSOLUTE ONE
% TO ZOOM IN/OUT WE ARE GOING TO MODIFY THE AXES LIMITS
stepx = (handles.axialxlim(1, 2) - handles.axialxlim(1, 1)) / 5;
stepy = (handles.axialylim(1, 2) - handles.axialylim(1, 1)) / 5;

axes(handles.axes3)
% zoom in
if zoom_value > handles.zoomvalue

    handles.axes3.XLim(1,1) = handles.axialxlim(1, 1) + (stepx / 2) * (zoom_value - 1);
    handles.axes3.XLim(1,2) = handles.axialxlim(1, 2) - (stepx / 2) * (zoom_value - 1);    
    handles.axes3.YLim(1,1) = handles.axialylim(1, 1) + (stepy / 2) * (zoom_value - 1);
    handles.axes3.YLim(1,2) = handles.axialylim(1, 2) - (stepy / 2) * (zoom_value - 1);
   
%     zoom(abs(handles.zoomvalue - zoom_value) + 1); % zoom(1) is not zooming
%     
%     % updating line thickness in the axes
%     axes(handles.axes3)
%     zoom out
%      zoom(zoom_value);
    
    axesHandlesToChildObjects = findobj(gca, 'Type', 'Line');
    if ~isempty(axesHandlesToChildObjects)
        for indexline = 1 : size(axesHandlesToChildObjects, 1)
            if strcmp(axesHandlesToChildObjects(indexline).LineStyle, 'none')
                set(axesHandlesToChildObjects, 'MarkerSize', round(zoom_value) * 10);
            else
                set(axesHandlesToChildObjects, 'linewidth', round(zoom_value));
            end
        end
    end
% end
end

% zoom out
if zoom_value < handles.zoomvalue
    handles.axes3.XLim(1,1) = handles.axialxlim(1, 1) + (stepx / 2) * (zoom_value - 1);
    handles.axes3.XLim(1,2) = handles.axialxlim(1, 2) - (stepx / 2) * (zoom_value - 1);
    handles.axes3.YLim(1,1) = handles.axialylim(1, 1) + (stepy / 2) * (zoom_value - 1);
    handles.axes3.YLim(1,2) = handles.axialylim(1, 2) - (stepy / 2) * (zoom_value - 1);
    %     oldzoomvalue01 = (handles.zoomvalue - 1) / (5 - 1);
    %
    %     newzoomvalue01 = (zoom_value - 1) / (5 - 1);
    %
    %     if newzoomvalue01 == 0
    %         zoom out
    %     else
    %         zoom(1 - abs(oldzoomvalue01 - newzoomvalue01));
    %     end
    %      zoom out
    %      zoom(zoom_value);
    %
    %     % updating line thickness in the axes
    %     axes(handles.axes3)
    %
    axesHandlesToChildObjects = findobj(gca, 'Type', 'Line');
    if ~isempty(axesHandlesToChildObjects)
        for indexline = 1 : size(axesHandlesToChildObjects, 1)
            if strcmp(axesHandlesToChildObjects(indexline).LineStyle, 'none')
                set(axesHandlesToChildObjects, 'MarkerSize', round(zoom_value) * 10);
            else
                set(axesHandlesToChildObjects, 'linewidth', round(zoom_value));
            end
        end
    end
end


% 
handles.zoomvalue = zoom_value;

guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function Slider_zoom_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Slider_zoom (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in refvolume_checkbox.
function refvolume_checkbox_Callback(hObject, eventdata, handles)
% hObject    handle to refvolume_checkbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of refvolume_checkbox

index_selected = get(handles.listbox1,'Value');

FIELDRTplotContours(handles, index_selected)


% --- Executes on slider movement.
function slider4_Callback(hObject, eventdata, handles)
% hObject    handle to slider4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end
