function [path_GS] = FIELDRTgetGSData(FIELDRT_path)
% FIELDRTgetGSData
% Imports the GS data into GUI listbox.
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

% Get the path of the directory to be selected for import.

if ispc
    path_GS = dir([FIELDRT_path '\GS\']);
    % Remove noises in the path_GS variable
    %   idx = ~cellfun('isempty',strfind({path_GS.name},'Oesophagus Case'));
    %   [i, j] = find(idx == 0);
    idx = cellfun('length', {path_GS.name});
    [i, j] = find(idx < 4);
    path_GS(j) = [];
    
    % To remove .DS_Store
    idx = find(ismember({path_GS.name}, '.DS_Store'));
    if idx ==0
        
    else
        path_GS(idx) = [];
    end
    
elseif isunix
    path_GS = dir([FIELDRT_path '/GS/']);
    % Remove noises in the path_GS variable
    %     idx = ~cellfun('isempty',strfind({path_GS.name},'Oesophagus Case'));
    %     [i, j] = find(idx == 0);
    idx = cellfun('length', {path_GS.name});
    [i, j] = find(idx < 5);
    path_GS(j) = [];
    
    % To remove .DS_Store
    idx = find(ismember({path_GS.name}, '.DS_Store'));
    if idx ==0
        
    else
        path_GS(idx) = [];
    end
    
end
end


