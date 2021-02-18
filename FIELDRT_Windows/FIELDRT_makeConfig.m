function FIELDRT_h = FIELDRT_makeConfig
%
% Retrieve FIELDRT configuration file (if this exists).
% Make FIELDRT configuration file (if this does not exists).
% FIELDRT configuration file is stored by default in $userHome/.FIELDRT/FIELDRT_Config.mat

% Emiliano Spezi 2013
% Edited by Concetta Piazzese Mar 2018
% import java.nio.file.Paths;
% import java.nio.file.Paths

% Retrieve basic environment data
userHome = char(java.lang.System.getProperty('user.home'));
FIELDRTDir = fullfile(userHome,'.FIELDRT');
FIELDRTConfigFilename = fullfile(FIELDRTDir, 'FIELDRT_Config.mat');

if exist(FIELDRTConfigFilename,'file')
    % Environment file exists 
    load(FIELDRTConfigFilename);
else
    % Environment file does not exist
    button = questdlg('FIELDRT config file does not exist? Do you want to create it.','FIELDRT environment','Yes','No','Yes');
    if strcmp(button,'Yes')
        FIELDRT_h.env.osName = char(java.lang.System.getProperty('os.name'));               %Operating system name
        FIELDRT_h.env.osArch = char(java.lang.System.getProperty('os.arch'));               %Operating system architecture
        FIELDRT_h.env.osVersion  = char(java.lang.System.getProperty('os.version'));        %Operating system version
        FIELDRT_h.env.fileSeparator = char(java.lang.System.getProperty('file.separator')); %File separator ("/" on UNIX)
%         FIELDRT_h.env.pathSeparator = char(java.lang.System.getProperty('path.separator')); %Path separator (":" on UNIX)
%         FIELDRT_h.env.lineSeparator = char(java.lang.System.getProperty('line.separator')); %Line separator ("\n" on UNIX)
        FIELDRT_h.env.userName = char(java.lang.System.getProperty('user.name'));           %User's account name
        FIELDRT_h.env.userHome = char(java.lang.System.getProperty('user.home'));           %User's home directory
        FIELDRT_h.env.userDir = char(java.lang.System.getProperty('user.dir'));             %User's current working directory
        FIELDRT_h.FIELDRTDir = FIELDRTDir;
        FIELDRT_h.FIELDRTConfigFilename = FIELDRTConfigFilename;
        if ~exist(FIELDRTDir,'dir')
            mkdir(FIELDRTDir);
        end
        save(FIELDRTConfigFilename,'FIELDRT_h');
    else
        FIELDRT_h = [];
    end
end



