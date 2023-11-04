%% ========== Setup ========== %% 
fclose('all'); clearvars; close all; clc;
addpath(genpath(pwd));

%% ========== Load LiDAR Data - VLP16 ========== %%
% ===== Select LiDAR File ===== %
lidarFile = '2023-06-09-15-51-42_Velodyne-VLP-16-Data.pcap';
lidarPath = 'C:\Users\user\Documents\NCKU-Data\competition\SG\data\';
%[lidarFile, lidarPath, ~] = uigetfile('*.pcap', 'Please select your LiDAR file');
if lidarFile == 0
    return;
end

[LiDAR_CFG, veloReader] = create_lidar_config(lidarPath, lidarFile);
t_lidar_start = PCAP_UTC_TO_GPS_TIME(seconds(veloReader.StartTime)*1000000, LiDAR_CFG);
t_lidar_end = PCAP_UTC_TO_GPS_TIME(seconds(veloReader.EndTime)*1000000, LiDAR_CFG);
dt = PCAP_UTC_TO_GPS_TIME(seconds(veloReader.Timestamps(2))*1000000, LiDAR_CFG) - ...
     PCAP_UTC_TO_GPS_TIME(seconds(veloReader.Timestamps(1))*1000000, LiDAR_CFG);
lidarRate = round(1/dt); % frequency (Hz)

if isequal(veloReader.DeviceModel, 'HDL64E')
    % Define vertical and horizontal resolution.
    LiDAR_CFG.vResolution = 64;
    LiDAR_CFG.hResolution = 1024;

    % Define vertical and horizontal field-of-view.
    LiDAR_CFG.vFoVUp = 2;     
    LiDAR_CFG.vFoVDown = -24.9; 
    LiDAR_CFG.vFoV = [LiDAR_CFG.vFoVUp, LiDAR_CFG.vFoVDown];
    LiDAR_CFG.hFoV = 360;

elseif isequal(veloReader.DeviceModel, 'VLP16')
    % Define vertical and horizontal resolution.
    LiDAR_CFG.vResolution = 16;
    LiDAR_CFG.hResolution = 1024;

    % Define vertical and horizontal field-of-view.
    LiDAR_CFG.vFoVUp = 15;     
    LiDAR_CFG.vFoVDown = -15; 
    LiDAR_CFG.vFoV = [LiDAR_CFG.vFoVUp, LiDAR_CFG.vFoVDown];
    LiDAR_CFG.hFoV = 360;

end

% ===== Specify Voxel Size for NDT ===== %
LiDAR_CFG.VoxelSizeNDT = 1;

%% Initial LiDAR Mounting Parameters - W.R.T. Right-Handed System
% ========== ENU-System ========== %
%     Z(Up)  Y(FWD)
%       |    /
%       |   /
%       |  /
%       | / 
%       *__ __ __ __ X(RIGHT)
%   (IMU's center)
% ================================ %

% ===== Lever-Arms [la_right; la_fwd; la_up] ===== %
r_lb = [0.0907; 1.1788; 1.2309];

% ===== Bore-Sight Angles [bs_right; bs_fwd; bs_up] ===== %
bs = [-1.0545; 0.3441518; -1.9687];
bsRotm = eul2rotmENU(bs);
C_bl = bsRotm.matrix;

%% ========== Load Integrated Navigation Data - TC-INS/GNSS ========== %%
initPOS = [23.0028474266, 120.214659453];
initPOS4WGS84 = [initPOS(1), initPOS(2), 38.070];
FileNameGT = 'REF_1.txt';
PathNameGT = 'C:\Users\user\Documents\NCKU-Data\competition\SG\data\';
%[FileNameGT, PathNameGT, ~] = uigetfile('*.txt', 'Please select your reference file (IE)');
f_INS = load([PathNameGT, FileNameGT]);
t_INS_start = f_INS(1, 1);
LLF_NED_INS = wgs2llf(f_INS, initPOS4WGS84);

disp(['Local Frame = ', sprintf('%7.4f, %7.4f', LLF_NED_INS(80000, 2), LLF_NED_INS(80000, 3))])
for i = 1:size(LLF_NED_INS, 1)
    if (LLF_NED_INS(i, 10) < 0)
        LLF_NED_INS(i, 10) = LLF_NED_INS(i, 10) + 360;
    end
    LLF_NED_INS(i, 2) = LLF_NED_INS(i, 2) + 2544814.048;
    LLF_NED_INS(i, 3) = LLF_NED_INS(i, 3) + 169492.873;
end
disp(['TWD97 Frame = ', sprintf('%10.4f, %11.4f', LLF_NED_INS(80000, 2), LLF_NED_INS(80000, 3))])

clearvars i f_ref;
%% ========== Read the First Point Cloud and Display It at the MATLAB ========== %%
ptCloud = readFrame(veloReader, 1);
disp(ptCloud)


%% ========== Visualize the Point Cloud ========== %% 
% ===== Create a Streaming Point Cloud Display Object ===== %
xlimits = [-45, 45];
ylimits = [-45, 45];
zlimits = [-10, 20];
lidarPlayer = pcplayer(xlimits, ylimits, zlimits);

% ===== Customize Player Axes Labels ===== %
xlabel(lidarPlayer.Axes, 'X (m)')
ylabel(lidarPlayer.Axes, 'Y (m)')
zlabel(lidarPlayer.Axes, 'Z (m)')
title(lidarPlayer.Axes, 'Lidar Sensor Data')

% ===== Visualize Point Cloud ===== %

for n = 1: 2
    ptCloud = readFrame(veloReader, n);
    view(lidarPlayer, ptCloud);

    pause(0.01)
end

%% ========== Build a Map Using Odometry ========== %%
% ===== Set Random Seed to Ensure Reproducibility ===== %
rng(0);

% ===== Create an Empty View Set ===== %
vSet = pcviewset;

% ===== Create a Figure for View Set Display ===== %
hFigBefore = figure('Name', 'View Set Display');
hAxBefore = axes(hFigBefore);

% ===== Initialize Point Cloud Processing Parameters ===== %
downsamplePercent = 0.1;
regGridSize = 1;

% ===== Initialize Transformations ===== %
absTform = rigidtform3d;  % Absolute transformation to reference frame
relTform = rigidtform3d;  % Relative transformation between successive scans

viewId = 1;
skipFrames = 10;
numFrames = veloReader.NumberOfFrames;
displayRate = 10; % Update display every 100 frames
start_ins_row = 1;
firstFrame = true;
currScan = 0;
t_INS = t_INS_start;

for n = 1: skipFrames: numFrames
    t_lidar = PCAP_UTC_TO_GPS_TIME(seconds(veloReader.Timestamps(n))*1000000, LiDAR_CFG);

    % Check an existing time of current LiDAR
    if (t_lidar < t_INS_start)
        continue;
    end

    % Time Synchronization:
    while (t_INS <= t_lidar)
        start_ins_row = start_ins_row + 1;
        t_INS = LLF_NED_INS(start_ins_row, 1);        
    end
    currScan = currScan + 1;

    % Interpolation
    interLocalPose = interpolatePOS(t_lidar, LLF_NED_INS, start_ins_row);
    tmpInterLocalpose = TF_NED2ENU(interLocalPose(:, 1: 10));
    LLF_ENU_INS.INT_SOL3(currScan, 1: 10) = tmpInterLocalpose;
    clearvars tmpInterLocalpose interLocalPose;
    nav_INS = getNavSolENU(LLF_ENU_INS.INT_SOL3, currScan);
    r_bn_INS = nav_INS.tran; % E, N, U
    R_n2b_INS = nav_INS.rotm;  % C_nb (n-frame: ENU to b-frame: Right-Front-Up)
    R_n2l_INS = C_bl * R_n2b_INS; % C_nl (n-frame: ENU to lidar-frame: Front-Left-Up)
    R_l2n_INS = R_n2l_INS'; % C_ln (lidar-frame: Front-Left-Up n-frame: ENU)
    r_ln_INS = (R_n2b_INS')*r_lb + r_bn_INS'; % Lidar position (from TC solution) in n-frame

    % Read point cloud
    ptCloudOrig = readFrame(veloReader,n);

    % Process point cloud
    ptCloud = helperProcessPointCloud(ptCloudOrig);

    % Visualize point cloud
    view(lidarPlayer, ptCloud);
    
    % Downsample the processed point cloud
    % ptCloud = pcdownsample(ptCloud, "random", downsamplePercent);
    
    if firstFrame
        absTform = rigidtform3d(R_l2n_INS, r_ln_INS');
        % Add first point cloud scan as a view to the view set
        % absTform = rigidtform3d(R_l2n_INS, r_ln_INS');
        vSet = addView(vSet, viewId, absTform, "PointCloud", ptCloudOrig);
        viewId = viewId + 1;
    
        % ptCloudPrev = pctransform(ptCloud, absTform);
        ptCloudPrev = ptCloud;
        firstFrame = false;
        continue;
    end
    
    % Use INS to estimate an initial transformation for registration
    delta_r_n = LLF_ENU_INS.INT_SOL3(currScan, 2: 4)' - LLF_ENU_INS.INT_SOL3(currScan-1, 2: 4)';
    nav_INS_prev = getNavSolENU(LLF_ENU_INS.INT_SOL3, currScan-1);
    R_n2l_INS_prev = C_bl * nav_INS_prev.rotm;  % C_nl_prev = C_bl*C_nb
    delta_r = R_n2l_INS_prev * delta_r_n;  % under lidar frame
    R_prev_now = (C_bl*R_n2b_INS) * R_n2l_INS_prev';
    R_now_prev = R_prev_now';
    initTform = rigidtform3d(R_now_prev, delta_r);
    
    % Compute rigid transformation that registers current point cloud with previous point cloud
    relTform = initTform;
    
    % Update absolute transformation to reference frame (first point cloud)
    absTform = rigidtform3d(R_l2n_INS, r_ln_INS');

    % Add current point cloud scan as a view to the view set
    vSet = addView(vSet, viewId, absTform, "PointCloud", ptCloudOrig);

    % Add a connection from the previous view to the current view, representing the relative transformation between them
    vSet = addConnection(vSet, viewId-1, viewId, relTform);
    viewId = viewId + 1;
    ptCloudPrev = ptCloud;
    initTform   = relTform;
    
    if n > 1 && mod(n, displayRate) == 1
        plot(vSet, "Parent", hAxBefore);
        drawnow update
    end
end

% ===== Save Output PCD ===== %
ptClouds = vSet.Views.PointCloud;
absPoses = vSet.Views.AbsolutePose;
mapGridSize = 0.2;
ptCloudMap = pcalign(ptClouds, absPoses, mapGridSize);

%output_pcd(ptCloudMap)

%% ========== Function ========== %%
function output_pcd(ptCloudMap)
    % Assuming ptCloud is your point cloud object and 'intensity' is a vector of intensity values
    minI = min(double(ptCloudMap.Intensity));
    maxI = max(double(ptCloudMap.Intensity));
    
    % Normalize intensity to 0-1
    normI = (double(ptCloudMap.Intensity)-minI) / (maxI-minI);
    
    % Map intensity to an RGB color (here we use a grayscale colormap, but others can be used)
    colorMap = uint8(255*repmat(normI, 1, 3));
    
    % Assign color to point cloud
    ptCloudMap.Color = colorMap;
    
    % Write to PLY file
    pcwrite(ptCloudMap, 'Odometry_vlp16_map.pcd');
end

function ptCloud = helperProcessPointCloud(ptCloudIn, method)
    % helperProcessPointCloud Process pointCloud to remove ground and ego vehicle
    arguments
        ptCloudIn (1,1) pointCloud
        method string {mustBeMember(method, ["planefit","rangefloodfill"])} = "rangefloodfill"
    end  
    isOrganized = ~ismatrix(ptCloudIn.Location);
    
    if (method == "rangefloodfill" && isOrganized) 
        % Segment ground using floodfill on range image
        groundFixedIdx = segmentGroundFromLidarData(ptCloudIn, "ElevationAngleDelta", 11);
    else
        % Segment ground as the dominant plane with reference normal vector pointing in positive z-direction
        maxDistance = 0.4;
        maxAngularDistance = 5;
        referenceVector = [0, 0, 1];
    
        [~, groundFixedIdx] = pcfitplane(ptCloudIn, maxDistance, referenceVector, maxAngularDistance);
    end
    
    if isOrganized
        groundFixed = false(size(ptCloudIn.Location, 1), size(ptCloudIn.Location, 2));
    else
        groundFixed = false(ptCloudIn.Count, 1);
    end
    groundFixed(groundFixedIdx) = true;
    
    % Segment ego vehicle as points within a given radius of sensor
    sensorLocation = [0, 0, 0];
    radius = 4;
    egoFixedIdx = findNeighborsInRadius(ptCloudIn, sensorLocation, radius);
    
    if isOrganized
        egoFixed = false(size(ptCloudIn.Location, 1),size(ptCloudIn.Location, 2));
    else
        egoFixed = false(ptCloudIn.Count, 1);
    end
    egoFixed(egoFixedIdx) = true;
    
    % Retain subset of point cloud without ground and ego vehicle
    if isOrganized
        indices = ~groundFixed & ~egoFixed;
    else
        indices = find(~groundFixed & ~egoFixed);
    end
    ptCloud = select(ptCloudIn, indices);

end

function helperMakeFigurePublishFriendly(hFig)
    if ~isempty(hFig) && isvalid(hFig)
        hFig.HandleVisibility = 'callback';
    end
end
