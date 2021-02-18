function [Overunder_imagetoshow, Jacc2D_cell] = FIELDRTstatsfiguresaving(varargin, range_total, struct_name_case, module)
% Creating and storing the figure to save in the pdf report

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
% Concetta Piazzese December 2019

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
handles.Case = varargin{1, 1}.Case;
handles.FIELDRTDataAttempts = varargin{1, 1}.FIELDRTDataAttempts;
handles.AttemptInformation = varargin{1, 1}.AttemptInformation;
handles.FIELDRTGSCases = varargin{1, 1}.FIELDRTGSCases;


%% Creating the 2D Jaccard distribution
f = figure('visible','off', 'Position', [1915 1915 900 300]);

% To have the figure with Minimal White Space
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1)*2.2;
bottom = outerpos(2) + ti(2)*2.2;
ax_width = outerpos(3) - ti(1)*2  - ti(3);
ax_height = outerpos(4) - ti(2)*3 - ti(4);
ax.Position = [left bottom ax_width ax_height];


%% Creating the traffic light image from the 2D
% Retrieving 2D Jaccard information
Jacc2D_cell = (handles.jaccard2Doutput); % Cell copy of handles.jaccard2Doutput
Jacc2D = cell2mat(handles.jaccard2Doutput);

% Finding the minimum and the maximum slice for each structures
minrange = min(range_total(:, 1));
maxrange = max(range_total(:, 2));

% Creating an RBG image all white
trafficlight_image = uint8(ones([size(Jacc2D), 3]) + 254); %RGB image

% Creating a matrix so to flag the slices to show
Overunder_imagetoshow = uint8(zeros(size(Jacc2D))); 

for index1 = 1 : size(Jacc2D, 1)
    for index2 = 1 : size(Jacc2D, 2)
        if Jacc2D(index1, index2) < 0.50 && (index2 >= range_total(index1, 1) || index2 <= range_total(index1, 2))
            trafficlight_image(index1, index2, 1) = 255;
            trafficlight_image(index1, index2, 2) = 0;
            trafficlight_image(index1, index2, 3) = 0;
            Overunder_imagetoshow(index1, index2) = 1; % flagging the image
        end
        if Jacc2D(index1, index2) >= 0.50 && Jacc2D(index1, index2) < 0.75
            trafficlight_image(index1, index2, 1) = 255;
            trafficlight_image(index1, index2, 2) = 255;
            trafficlight_image(index1, index2, 3) = 0;
        end
        if Jacc2D(index1, index2) >= 0.75
            trafficlight_image(index1, index2, 1) = 0;
            trafficlight_image(index1, index2, 2) = 255;
            trafficlight_image(index1, index2, 3) = 0;
        end
        
        % for slices traced by the user but outside the GS range
        if Jacc2D(index1, index2) == 0 && (index2 < range_total(index1, 1) || index2 > range_total(index1, 2))
            trafficlight_image(index1, index2, 1) = 202;
            trafficlight_image(index1, index2, 2) = 204;
            trafficlight_image(index1, index2, 3) = 206;
            Overunder_imagetoshow(index1, index2) = 1; % flagging the image
        end
        
        if isnan(Jacc2D(index1, index2)) && (index2 < range_total(index1, 1) || index2 > range_total(index1, 2))
            trafficlight_image(index1, index2, 1) = 134;
            trafficlight_image(index1, index2, 2) = 136;
            trafficlight_image(index1, index2, 3) = 138;
            Overunder_imagetoshow(index1, index2) = 0; % flagging the image
            
            Jacc2D(index1, index2) = 0;
            Jacc2D_cell{index1, index2} = 0;
        end
    end
end

% 

trafficlight_image_toshow = trafficlight_image(:, minrange : maxrange, :);

%% Creating the traffic light image
f = figure('visible','off', 'Position', [1915 1915 1800 350]);

% To have the figure with Minimal White Space
ax = gca;
outerpos = ax.OuterPosition;
ti = ax.TightInset;
left = outerpos(1) + ti(1)*5;
bottom = outerpos(2) + ti(2)*4;
ax_width = outerpos(3) - ti(1)*6- ti(3);
ax_height = outerpos(4) - ti(2)*6.5 - ti(4);
ax.Position = [left bottom ax_width ax_height];


himshow = imshow(trafficlight_image_toshow);
axis ij % To have the 0 at the left bottom corner

% Setting some image properties
set(gca, 'TickLength', [0 0]) % to remove/edit the grid axis short-line

% Labels and position of the labels on y-axis
ynames = struct_name_case';
ytick = 1 : size(ynames, 1);

% Labels and position of the labels on x-axis
xtick = 1 : size(trafficlight_image_toshow, 2);
for indextick = 1 : size(trafficlight_image_toshow, 2)
    xnames(1, indextick) = {num2str(minrange + xtick(indextick) - 1)};
end

% setting the x-,y-lim
xlim([0.5 (size(trafficlight_image_toshow, 2) + 0.5)]);
ylim([0.5 (size(ynames, 1) + 0.5)])

set(gca,'FontSize', 22, 'xtick', xtick,'xticklabel', xnames, 'ytick', ytick,'yticklabel', ynames)

axis on

hold on

% Adding grid lines
M = size(trafficlight_image_toshow, 1);
N = size(trafficlight_image_toshow, 2);

for k = 0.5 : (M + 0.5)
    x = [-1 (N + 1)];
    y = [k k];
    plot(x,y,'Color','k','LineStyle','-');
    plot(x,y,'Color','k','LineStyle','-');
end

for k = 0.5 : (N + 0.5)
    x = [k k];
    y = [-1 (M + 1)];
    plot(x,y,'Color','k','LineStyle','-');
    plot(x,y,'Color','k','LineStyle','-');
end

caption = sprintf('2D Jaccard');
title(caption, 'FontSize', 40, 'Color' , 'k');
set(gcf,'color','w');

xlabel('Slice', 'FontSize', 28, 'FontWeight', 'bold'); % x-axis label
% ylabel ('2D Jaccard', 'FontSize', 28, 'FontWeight', 'bold'); % y-axis label



if ispc
    % saveas(f, [handles.FIELDRT.env.userDir '\FIELDRT_Reports\Images\' module '\'  struct_name_case{index1} '\' struct_name_temp int2str(index2)] , 'tiff')
    print(f, [handles.FIELDRT.env.userDir '\FIELDRT_Reports\Images\' module '\Stats\Traffic_light_stats'] , '-dtiffn')
else
    % saveas(f, [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/' struct_name_case{index1} '/' struct_name_temp int2str(index2)] , 'tiff')
    print(f, [handles.FIELDRT.env.userDir '/FIELDRT_Reports/Images/' module '/Stats/Traffic_light_stats'] , '-dtiffn')
end


% Closing the figure
clear f

end
