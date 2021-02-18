function FIELDRTquit_Callback(hObject, eventdata, handles)
% FIELDRTquitbutton_Callback
% Close FIELDRT viewer
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
% Concetta Piazzese Januart 2019

close FIELDRTResultvisualization

cd(handles.FIELDRT_h.env.userDir) % changes the current folder to newFolder.

% Reopening FIELDRT
FIELDRT

end
