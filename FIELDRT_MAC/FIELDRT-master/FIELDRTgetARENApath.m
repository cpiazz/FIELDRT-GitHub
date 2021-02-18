function [FIELDRT_path] = FIELDRTgetFIELDRTpath(str)
% FIELDRTgetFIELDRTpath
% Gets the FIELDRT folder path.
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
    indV = find(str == '\');
    ind = max(indV);
    FIELDRT_path = [str(1:ind) 'FIELDRT-master\'];
elseif isunix
    indV = find(str == '/');
    ind = max(indV);
    FIELDRT_path = [str(1:ind) 'FIELDRT-master/'];
end

% Adding the folder to path
addpath(genpath(FIELDRT_path));

end





