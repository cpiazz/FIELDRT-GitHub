function [planC] = FIELDRTrenameGSstructures(planC_temp)
% FIELDRTrenameGSstructures
% Renames the GS structures in the variable planC.
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

planC = planC_temp;
for k = 1 : size(planC_temp{1, 4}, 2)
    planC{1, 4}(k).structureName = [planC_temp{1, 4}(k).structureName '_GS'];
end
end

