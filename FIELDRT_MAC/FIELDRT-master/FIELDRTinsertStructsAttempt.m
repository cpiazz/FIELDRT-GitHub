function [planC, DateVector, AttemptInformation, status, Datatopass] = FIELDRTinsertStructsAttempt(planC, FIELDRT_path_FIELDRT, FIELDRTDataAttempts, handles, status)
% FIELDRTinsertStructsAttempt
% Insert the structures of a DICOM file in the variable planC.
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
% CPiazzese

% import org.dcm4che2.io.*
% import org.dcm4che2.data.*

warning off

%Start diary and timer for import log.
startTime = now;
tmpFileName = tempname;
diary(tmpFileName);

t = [datetime('now')];
DateVector = datevec(t);

%Load study containing the structure(s) to be inserted.
[fileName, dirName]= uigetfile('*.dcm', ...
    'Select the study containing the structure(s) to insert.');  %at position 100, 100 in pixels.
dcmFileName = fullfile(dirName,fileName);

if fileName == 0 % the user press cancel
    DateVector = [];
    AttemptInformation = [];
    Datatopass = [];
    status = 4;
    return
end

AttemptInformation.FileName = fileName;
AttemptInformation.FiledirName = dirName;

% % To avoid the user to close the window
% set(handles.figure1,'CloseRequestFcn', ' ');

AttemptInformation.numAttempt = 0;
dataAttempt = [];
Datatopass = [];

%% Checking if the attempt has been already loaded (now using the SOP Instance UID)
% for index = 1 : size(FIELDRTDataAttempts, 2)
%     if (strcmp(fileName, FIELDRTDataAttempts(index).AttemptName) && strcmp(dirName, FIELDRTDataAttempts(index).AttemptPath))
%         dataAttempt = [dataAttempt; index];
%         AttemptInformation.numAttempt = AttemptInformation.numAttempt + 1;
%         AttemptInformation.index = index; % to save the index of the attempt in the FIELDRTDataAttempts variable
%     end
% end
% Loading the DICOM info of the attempt and extracting the SOPInstanceUID tag
infodicom = dicominfo([AttemptInformation.FiledirName AttemptInformation.FileName]);

for index = 1 : size(FIELDRTDataAttempts, 2)
    if strcmp(infodicom.SOPInstanceUID, FIELDRTDataAttempts(index).SOPInstanceUID)
        AttemptInformation.numAttempt = AttemptInformation.numAttempt + 1;
        AttemptInformation.index = index; % to save the index of the attempt in the FIELDRTDataAttempts variable
    end
end


% Creating the report folder
if ispc
    if ~exist([handles.FIELDRT_h.env.userDir '\FIELDRT_Reports'] ,'dir')
        mkdir([handles.FIELDRT_h.env.userDir '\FIELDRT_Reports']);
    end
else
    if ~exist([handles.FIELDRT_h.env.userDir '/FIELDRT_Reports'] ,'dir')
        mkdir([handles.FIELDRT_h.env.userDir '/FIELDRT_Reports']);
    end
end

% Coding for the status
% status = 1 first attempt correctly loaded
% status = 2 attempt already submitted



if AttemptInformation.numAttempt == 0
    % Attempt has never been analysed
    h = msgbox('This is your first attempt');
    waitfor(h);
    
    status = 1;
    
    AttemptInformation.SOPInstanceUID = infodicom.SOPInstanceUID; % saving the SOP Instance UID of the attempt
else
    % Attempt has been analyzed
    %     choice = questdlg(['This attempt has been submitted ' int2str(AttemptInformation.numAttempt) ' times. Click to submit new attempt'], ...
    %         ' ', 'Submit new attempt', ...
    %         'View results of this attempt', 'View results of this attempt');
    choice = questdlg(['This attempt has been submitted'], ...
        ' ', 'Submit new attempt', ...
        'View results of this attempt', 'View results of this attempt');
    % Handle response
    switch choice
        case 'Submit new attempt'
            [fileName, dirName]= uigetfile('*.dcm', ...
                'Select the study containing the structure(s) to insert.');  %at position 100, 100 in pixels.
            dcmFileName = fullfile(dirName,fileName);
            if fileName == 0
                status = 4; % FIELDRT will be restarted
                return
            else
                
                AttemptInformation.FileName = fileName;
                AttemptInformation.FiledirName = dirName;
                
                AttemptInformation.numAttempt = 0;
                dataAttempt = [];
                
                %% Checking if the attempt has been already loaded (now using the SOP Instance UID)
                infodicom = dicominfo([AttemptInformation.FiledirName AttemptInformation.FileName]);
                
                for index = 1 : size(FIELDRTDataAttempts, 2)
                    if strcmp(infodicom.SOPInstanceUID, FIELDRTDataAttempts(index).SOPInstanceUID)
                        AttemptInformation.numAttempt = AttemptInformation.numAttempt + 1;
                        AttemptInformation.index = index; % to save the index of the attempt in the FIELDRTDataAttempts variable
                    end
                end
                
                if AttemptInformation.numAttempt == 0
                    % Attempt has never been analysed
                    h = msgbox('This is your first attempt');
                    waitfor(h);
                    
                    AttemptInformation.SOPInstanceUID = infodicom.SOPInstanceUID; % saving the SOP Instance UID of the attempt
                    status = 1;
                else
                    h1 = msgbox('You have already submitted this attempt');
                    waitfor(h1);
                    %                     if ispc
                    %
                    %                         % Loading the saved attempt
                    %                         load([handles.FIELDRT_h.env.userDir '\FIELDRT-master\MatAttempts\' fileName(1 : end - 4) '_' ...
                    %                             int2str(FIELDRTDataAttempts(dataAttempt(1, 1)).Date.Day) '_' int2str(FIELDRTDataAttempts(dataAttempt(1, 1)).Date.Month) '_' ...
                    %                             int2str(FIELDRTDataAttempts(dataAttempt(1, 1)).Date.Year) '_' int2str(FIELDRTDataAttempts(dataAttempt(1, 1)).Time.Hour) '_' ...
                    %                             int2str(FIELDRTDataAttempts(dataAttempt(1, 1)).Time.Minutes) '.mat']);
                    %
                    %                         %                     % Opening the pdf report and the result viewer
                    %                         %                     winopen([handles.FIELDRT_h.env.userDir '/Reports']);
                    %
                    %                         %                     set(handles.figure1,'CloseRequestFcn', 'closereq;')% To avoid quitbutton_Callback
                    %                         %
                    %                         %                     waitfor(FIELDRTResultvisualization(Datatopass));
                    %                         %
                    %                         %                     cd(handles.FIELDRT_h.env.userDir) % changes the current folder to newFolder.
                    %                         %                     FIELDRT
                    %                         status = 2;
                    %
                    %                     elseif isunix
                    %                         % Loading the saved attempt
                    %                         load([handles.FIELDRT_h.env.userDir '/FIELDRT-master/MatAttempts/' fileName(1 : end - 4) '_' ...
                    %                             int2str(FIELDRTDataAttempts(dataAttempt(1, 1)).Date.Day) '_' int2str(FIELDRTDataAttempts(dataAttempt(1, 1)).Date.Month) '_' ...
                    %                             int2str(FIELDRTDataAttempts(dataAttempt(1, 1)).Date.Year) '_' int2str(FIELDRTDataAttempts(dataAttempt(1, 1)).Time.Hour) '_' ...
                    %                             int2str(FIELDRTDataAttempts(dataAttempt(1, 1)).Time.Minutes) '.mat']);
                    %
                    %                         %                     % Opening the pdf report and the result viewer
                    %                         %                     inp = [handles.FIELDRT_h.env.userDir '\Reports'];
                    %                         %                     syscmd = ['open "', inp, '" &'];
                    %                         %                     %                     disp(['Running the following in the Terminal: "', syscmd,'"']);
                    %                         %                     system(syscmd);
                    %
                    %                         %                     set(handles.figure1,'CloseRequestFcn', 'closereq;')% To avoid quitbutton_Callback
                    %                         %
                    %                         %                     waitfor(FIELDRTResultvisualization(Datatopass));
                    %                         %
                    %                         %                     cd(handles.FIELDRT_h.env.userDir) % changes the current folder to newFolder.
                    %                         %                     FIELDRT
                    %                         status = 2;
                    %                     end
                    status = 4;
                    return
                end
            end
            
        case 'View results of this attempt'
            if ispc
                % Loading the saved attempt
                load([handles.FIELDRT_h.env.userDir '\FIELDRT-master\MatAttempts\' FIELDRTDataAttempts(index).AttemptName(1 : end - 4) '_' ...
                    int2str(FIELDRTDataAttempts(index).Date.Day) '_' int2str(FIELDRTDataAttempts(index).Date.Month) '_' ...
                    int2str(FIELDRTDataAttempts(index).Date.Year) '_' int2str(FIELDRTDataAttempts(index).Time.Hour) '_' ...
                    int2str(FIELDRTDataAttempts(index).Time.Minutes) '.mat']);
                
                %                 winopen([handles.FIELDRT_h.env.userDir '/Reports']);
                %                 load(FIELDRT_path_FIELDRT, '-mat');
                %                 % Opening the pdf report and the result viewer
                %                 inp = [handles.FIELDRT_h.env.userDir '\Reports'];
                %                 syscmd = ['open "', inp, '" &'];
                %                 %                     disp(['Running the following in the Terminal: "', syscmd,'"']);
                %                 system(syscmd);
                
                %                 set(handles.figure1,'CloseRequestFcn', 'closereq;')% To avoid quitbutton_Callback
                %
                %                 waitfor(FIELDRTResultvisualization(Datatopass));
                %
                %                 cd(handles.FIELDRT_h.env.userDir) % changes the current folder to newFolder.
                %                 FIELDRT
                status = 2;
                
            elseif isunix
                % Loading the saved attempt
                load([handles.FIELDRT_h.env.userDir '/FIELDRT-master/MatAttempts/' FIELDRTDataAttempts(index).AttemptName(1 : end - 4) '_' ...
                    int2str(FIELDRTDataAttempts(index).Date.Day) '_' int2str(FIELDRTDataAttempts(index).Date.Month) '_' ...
                    int2str(FIELDRTDataAttempts(index).Date.Year) '_' int2str(FIELDRTDataAttempts(index).Time.Hour) '_' ...
                    int2str(FIELDRTDataAttempts(index).Time.Minutes) '.mat']);
                
                %                 inp = [handles.FIELDRT_h.env.userDir '\Reports'];
                %                 syscmd = ['open "', inp, '" &'];
                %                 %                     disp(['Running the following in the Terminal: "', syscmd,'"']);
                %                 system(syscmd);
                
                %                 % Opening the pdf report and the result viewer
                %                 inp = [handles.FIELDRT_h.env.userDir '\Reports'];
                %                 syscmd = ['open "', inp, '" &'];
                %                     disp(['Running the following in the Terminal: "', syscmd,'"']);
                %                 system(syscmd);
                %
                %                 set(handles.figure1,'CloseRequestFcn', 'closereq;')% To avoid quitbutton_Callback
                %
                %                 waitfor(FIELDRTResultvisualization(Datatopass));
                %
                %                 cd(handles.FIELDRT_h.env.userDir) % changes the current folder to newFolder.
                %                 FIELDRT
                
                status = 2;
                
                
            end
            
            return
    end
    % To uncomment before releasing the software
    if isempty(choice)
        status = 4;
        return
    end
end

%% Start loading
% Deactivation of some buttons
handles.statictext.String = sprintf(['Loading attempt...\n Please wait']);
handles.loadbutton.Enable = 'off';
handles.listbox.Enable = 'off';
handles.loadbutton.Enable = 'off';
handles.analyseandresultsbutton.Enable = 'off';
handles.aboutbutton.Enable = 'off';
handles.helpbutton.Enable = 'off';
handles.quitbutton.Enable = 'off';

% Sets env variables necessary for operation of ML_DICOM.
initFlag = init_ML_DICOM;


set(handles.figure1,'CloseRequestFcn', @closereq) % To avoid quitbutton_Callback

indexS = planC{end};

[dcmObj, isDcm] = scanfile_mldcm(dcmFileName);
dcmdirS = [];
if isDcm
    dcmdirS = dcmdir_add(dcmFileName, dcmObj, dcmdirS);
    dcmObj.clear;
else
    error('Invalid DICOM data')
end

dataS = populate_planC_field('structures', dcmdirS.PATIENT);

% Tolerance to determine oblique scan (think about passing it as a
% parameter in future)
numScans = length(planC{indexS.scan});
obliqTol = 1e-3;
isObliqScanV = ones(1,numScans);

% Check scan to associate the strucutres
scanUIDc = {planC{indexS.scan}.scanUID};
scanTypesC = {};
for i = 1 : numScans
    
    scanTypesC{i} = [num2str(i) '.  ' planC{indexS.scan}(i).scanType];
    
    if isfield(planC{indexS.scan}(i).scanInfo(1),'DICOMHeaders') && ...
            ~isempty(planC{indexS.scan}(i).scanInfo(1).DICOMHeaders)
        ImageOrientationPatientV = planC{indexS.scan}(i).scanInfo(1).DICOMHeaders.ImageOrientationPatient;
    else
        ImageOrientationPatientV = [];
    end
    % Check for obliqueness
    if ~isempty(ImageOrientationPatientV) && max(abs((abs(ImageOrientationPatientV) - [1 0 0 0 1 0]'))) <= obliqTol
        isObliqScanV(i) = 0;
    end
    
end

% Return if not a valid RTStruct file
if isempty(dataS)
    CERRStatusString('File not valid.')
    return
end

% So that all the structures have the same associatedScan.
for indexstruct = 1 : size(dataS, 2)
    dataS(indexstruct).assocScanUID = planC{1, 4}(1).assocScanUID;
end

if ~ismember(dataS(1).assocScanUID,scanUIDc)
    %     scanInd = listdlg('PromptString','Select Scan to associate structures',...
    %         'SelectionMode','single',...
    %         'ListString',scanTypesC);
    scanInd = 1;
    if ~isempty(scanInd)
        [dataS.assocScanUID] = deal(scanUIDc{scanInd});
    else
        return
    end
end
numStructs = length(planC{indexS.structures});
for i=1:length(dataS)
    dataS(i) = sortStructures(dataS(i), isObliqScanV, planC);
    %     colorNum = numStructs + i;
    %     if isempty(dataS(i).structureColor)
    %         color = stateS.optS.colorOrder( mod(colorNum-1, size(stateS.optS.colorOrder,1))+1,:);
    %         dataS(i).structureColor = color;
    %     end
end
if ~isempty(planC{indexS.structures})
    for strNum = 1:length(dataS)
        planC{indexS.structures} = dissimilarInsert(planC{indexS.structures},dataS(strNum),length(planC{indexS.structures})+1);
    end
else
    planC{indexS.structures} = dataS;
end

editStructNumV = (numStructs+1):length(planC{indexS.structures});


planC = getRasterSegs(planC, 6);
%planC = setUniformizedData(planC);
for strNum = editStructNumV
    planC = updateStructureMatrices(planC, strNum);
end

% if isfield(stateS, 'CERRFile')
%     stateS.structsChanged = 1;
%     CERRRefresh
% end

% CERRStatusString('Done inserting structures.')

clear temp_planC;
