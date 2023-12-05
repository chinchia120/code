%% =============== Function of helperVisualizeEgoView =============== %%
function player = helperVisualizeEgoView(ptCloud)
    % Create a pcplayer object
    xlimits = ptCloud.XLimits;
    ylimits = ptCloud.YLimits;
    zlimits = ptCloud.ZLimits;
    player = pcplayer(xlimits, ylimits, zlimits);

    % Turn off axes lines
    axis(player.Axes, 'off');

    % Set up camera to show ego view
    camproj(player.Axes, 'perspective');
    camva(player.Axes, 0);
    campos(player.Axes, [0, 0, 0]);
    camtarget(player.Axes, [-1, 0, 0]);
    % view(player,ptCloud);

    %Set up a transformation to rotate by 5 degrees
    theta = 5;

    eulerAngles = [0, 0, theta];
    translation = [0, 0, 0];
    rotateByTheta = rigidtform3d(eulerAngles, translation);

    for n = 0 : theta : 359
        % Rotate point cloud by theta
        ptCloud = pctransform(ptCloud, rotateByTheta);

        % Display point cloud
        view(player, ptCloud);

        pause(0.05)
    end
end