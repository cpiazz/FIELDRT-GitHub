function FIELDRTMacopen(inp)
% Opens a file or directory using the OPEN terminal utility on the MAC.
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
% Concetta Piazzese March 2020

% if strcmpi('.',inp)
%     inp = pwd;
% end

syscmd = ['open "', inp, '"'];
% disp(['Running the following in the Terminal: "', syscmd,'"']);
system(syscmd);