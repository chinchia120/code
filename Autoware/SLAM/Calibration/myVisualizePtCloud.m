%==========================================================================
% Function: Visualize point cloud both before and after registration
% Developed & Coded by: Surachet Srinara, PhD candidate
% Date: September 5, 2021 - 02:28PM
% Copyright@Chet2021
%==========================================================================
function myVisualizePtCloud(currScanFrame, fixedScan, movingScanBefore, movingScanAfter, colorBefore, colorAfter)
% Change color for fixed scan frame:
color_ = ones(fixedScan.Count,3)*200;
fixedScan.Color = uint8(color_);

% Colorize point cloud:
if isequal(colorBefore,'r')
    pcColorBefore = [1 0 0];
elseif isequal(colorBefore,'g')
    pcColorBefore = [0 1 0];
elseif isequal(colorBefore,'b')
    pcColorBefore = [0 0 1];
else
    pcColorBefore = [1 1 0];
end

if isequal(colorAfter,'r')
    pcColorAfter = [1 0 0];
elseif isequal(colorAfter,'g')
    pcColorAfter = [0 1 0];
elseif isequal(colorAfter,'b')
    pcColorAfter = [0 0 1];
else
    pcColorAfter = [1 1 0];
end

% Change color for (before) moving scan frame:
thisPtCloud = pointCloud([movingScanBefore.Location(:,1),movingScanBefore.Location(:,2),movingScanBefore.Location(:,3)]);
cmatrix3 = ones(size(thisPtCloud.Location)).*pcColorBefore;
thisPtCloud = pointCloud([movingScanBefore.Location(:,1),movingScanBefore.Location(:,2),movingScanBefore.Location(:,3)], ...
                           'Color',cmatrix3);
movingScanBefore = thisPtCloud;  
clearvars thisPtCloud

% Change color for (after) moving scan frame:
thisPtCloud = pointCloud([movingScanAfter.Location(:,1),movingScanAfter.Location(:,2),movingScanAfter.Location(:,3)]);
cmatrix3 = ones(size(thisPtCloud.Location)).*pcColorAfter;
thisPtCloud = pointCloud([movingScanAfter.Location(:,1),movingScanAfter.Location(:,2),movingScanAfter.Location(:,3)], ...
                           'Color',cmatrix3);
movingScanAfter = thisPtCloud;  
clearvars thisPtCloud

% Point cloud show pair:
hFigAlign = figure;
axAlign1 = subplot(1, 2, 1,'Color', [0, 0, 0], 'Parent', hFigAlign);
pcshowpair(fixedScan, movingScanBefore, 'BlendFactor', 0, 'Parent', axAlign1);
title(axAlign1, sprintf('Before NDT-Registration [Scan No.: %0.0f]', currScanFrame))
% title(axAlign1, 'Before NDT-Registration');
xlabel('X (m.)');
ylabel('Y (m.)');
zlabel('Z (m.)');
view(axAlign1, 2);
axAlign2 = subplot(1, 2, 2,'Color', [0, 0, 0], 'Parent', hFigAlign);
pcshowpair(fixedScan, movingScanAfter, 'BlendFactor', 0, 'Parent', axAlign2);
title(axAlign2, sprintf('After NDT-Registration [Scan No.: %0.0f]', currScanFrame))
% title(axAlign2, 'After NDT-Registration');
xlabel('X (m.)');
ylabel('Y (m.)');
zlabel('Z (m.)');
view(axAlign2, 2);

end