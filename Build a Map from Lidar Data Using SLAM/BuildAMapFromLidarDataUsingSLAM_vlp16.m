%% Build a Map from Lidar Data Using SLAM
% This example shows how to process 3-D lidar data from a sensor mounted on 
% a vehicle to progressively build a map and estimate the trajectory of a vehicle 
% using simultaneous localization and mapping (SLAM). In addition to 3-D lidar 
% data, an inertial navigation sensor (INS) is also used to help build the map. 
% Maps built this way can facilitate path planning for vehicle navigation or can 
% be used for localization.
%% Overview
% The <docid:driving_ug#mw_f17b046b-90a0-40c1-a57f-d1a49ceb6ce0 Build a Map 
% from Lidar Data> example uses 3-D lidar data and IMU readings to progressively 
% build a map of the environment traversed by a vehicle. While this approach builds 
% a locally consistent map, it is suitable only for mapping small areas. Over 
% longer sequences, the drift accumulates into a significant error. To overcome 
% this limitation, this example recognizes previously visited places and tries 
% to correct for the accumulated drift using the graph SLAM approach.

fclose('all'); clearvars; close all; clc
addpath(genpath(pwd))
%% Load LiDAR Data:

% == First LiDAR: VLP-16 (Front-LiDAR)
[lidarFile, lidarPath, ~] = uigetfile('*.pcap', 'Please select your LiDAR file');
if lidarFile == 0
    return;
end
[LiDAR_CFG, veloReader] = create_lidar_config(lidarPath, lidarFile);
t_lidar_start = PCAP_UTC_TO_GPS_TIME(seconds(veloReader.StartTime)*1000000, LiDAR_CFG);
t_lidar_end = PCAP_UTC_TO_GPS_TIME(seconds(veloReader.EndTime)*1000000, LiDAR_CFG);
dt = (PCAP_UTC_TO_GPS_TIME(seconds(veloReader.Timestamps(2))*1000000, LiDAR_CFG)-PCAP_UTC_TO_GPS_TIME(seconds(veloReader.Timestamps(1))*1000000, LiDAR_CFG));
lidarRate = round(1/dt); % Hz
if isequal(veloReader.DeviceModel, 'HDL64E')
    % Define vertical and horizontal resolution.
    LiDAR_CFG.vResolution = 64;
    LiDAR_CFG.hResolution = 1024;
    % Define vertical and horizontal field-of-view.
    LiDAR_CFG.vFoVUp      = 2;     
    LiDAR_CFG.vFoVDown    = -24.9; 
    LiDAR_CFG.vFoV        = [LiDAR_CFG.vFoVUp LiDAR_CFG.vFoVDown];
    LiDAR_CFG.hFoV        = 360;
elseif isequal(veloReader.DeviceModel, 'VLP16')
    % Define vertical and horizontal resolution.
    LiDAR_CFG.vResolution = 16;
    LiDAR_CFG.hResolution = 1024;
    % Define vertical and horizontal field-of-view.
    LiDAR_CFG.vFoVUp      = 15;     
    LiDAR_CFG.vFoVDown    = -15; 
    LiDAR_CFG.vFoV        = [LiDAR_CFG.vFoVUp LiDAR_CFG.vFoVDown];
    LiDAR_CFG.hFoV        = 360;
end
% Specify voxel size for NDT:
LiDAR_CFG.VoxelSizeNDT = 1;
%% Initial LiDAR Mounting Parameters: W.R.T. right-handed system

% ********** ENU-System **********
%     Z(Up)  Y(FWD)
%       |    /
%       |   /
%       |  /
%       | / 
%       *__ __ __ __ X(RIGHT)
%   (IMU's center)
% *******************************
% == Lever-arms: [la_right; la_fwd; la_up]
r_lb   = [0.0907;1.1788;1.2309];

% == Bore-sight angles: [bs_right; bs_fwd; bs_up]
bs    = [-1.0545; 0.3441518; -1.9687];
bsRotm = eul2rotmENU(bs);
C_bl  = bsRotm.matrix;
%%  Load Integrated Navigation Data: TC-INS/GNSS:

initPOS                         = [23.0028474266, 120.214659453];
initPOS4WGS84                   = [initPOS(1), initPOS(2), 38.070];
[FileNameGT, PathNameGT, ~]     = uigetfile('*.txt', 'Please select your reference file (IE)');
f_INS                           = load([PathNameGT FileNameGT]);
t_INS_start                            = f_INS(1, 1);
LLF_NED_INS                     = wgs2llf(f_INS, initPOS4WGS84); % LLH2NED
for i = 1:size(LLF_NED_INS,1)
    if (LLF_NED_INS(i,10) < 0)
        LLF_NED_INS(i,10) = LLF_NED_INS(i,10) + 360;
    else
        continue;
    end
end
clearvars i f_ref;
%% 
% Read the first point cloud and display it at the MATLAB® command prompt

ptCloud = readFrame(veloReader,1);
disp(ptCloud)
%% 
% Visualize the point clouds using <docid:vision_ref#busm5az-7 |pcplayer|>, 
% a streaming point cloud display. The vehicle traverses a path consisting of 
% two loops. In the first loop, the vehicle makes a series of turns and returns 
% to the starting point. In the second loop, the vehicle makes a series of turns 
% along another route and again returns to the starting point.

% Specify limits of the player
xlimits = [-45 45]; % meters
ylimits = [-45 45];
zlimits = [-10 20];

% Create a streaming point cloud display object
lidarPlayer = pcplayer(xlimits, ylimits, zlimits);

% Customize player axes labels
xlabel(lidarPlayer.Axes, 'X (m)')
ylabel(lidarPlayer.Axes, 'Y (m)')
zlabel(lidarPlayer.Axes, 'Z (m)')

title(lidarPlayer.Axes, 'Lidar Sensor Data')
% Skip evey other frame since this is a long sequence
skipFrames  = 5;
% numFrames   = veloReader.NumberOfFrames;
numFrames   = 10;
for n = 1 : skipFrames : numFrames
    
    ptCloud = readFrame(veloReader,n);
    
    % Visualize point cloud
    view(lidarPlayer, ptCloud);
       
    pause(0.01)
end
%% Build a Map Using Odometry
% First, use the approach explained in the <docid:driving_ug#mw_f17b046b-90a0-40c1-a57f-d1a49ceb6ce0 
% Build a Map from Lidar Data> example to build a map. The approach consists of 
% the following steps:
%% 
% * *Align lidar scans:* Align successive lidar scans using a point cloud registration 
% technique. This example uses <docid:vision_ref#mw_004d7ab5-9956-4970-8919-107b4a2d2aa1 
% |pcregisterndt|> for registering scans. By successively composing these transformations, 
% each point cloud is transformed back to the reference frame of the first point 
% cloud.
% * *Combine aligned scans:* Generate a map by combining all the transformed 
% point clouds.
%% 
% This approach of incrementally building a map and estimating the trajectory 
% of the vehicle is called _odometry_.
% 
% Use a <docid:vision_ref#mw_c6968c3b-08dd-4aaa-94e7-b6c0e0185f8d |pcviewset|> 
% object to store and manage data across multiple views. A view set consists of 
% a set of connected views.
%% 
% * Each view stores information associated with a single view. This information 
% includes the absolute pose of the view, the point cloud sensor data captured 
% at that view, and a unique identifier for the view. Add views to the view set 
% using <docid:vision_ref#bu5a53u-1 |addView|>.
% * To establish a connection between views use <docid:vision_ref#mw_a70300e5-ecc6-4ecf-a6f8-9825494ef0f1 
% |addConnection|>. A connection stores information like the relative transformation 
% between the connecting views, the uncertainty involved in computing this measurement 
% (represented as an information matrix) and the associated view identifiers.
% * Use the <docid:vision_ref#mw_bad79ccd-2ea8-42e9-8e9f-4df662dec050 |plot|> 
% method to visualize the connections established by the view set. These connections 
% can be used to visualize the path traversed by the vehicle.

% hide(lidarPlayer)

% Set random seed to ensure reproducibility
rng(0);

% Create an empty view set
vSet = pcviewset;

% Create a figure for view set display
hFigBefore = figure('Name', 'View Set Display');
hAxBefore = axes(hFigBefore);

% Initialize point cloud processing parameters
downsamplePercent = 0.1;
regGridSize       = 1;

% Initialize transformations
absTform   = rigidtform3d;  % Absolute transformation to reference frame
relTform   = rigidtform3d;  % Relative transformation between successive scans

viewId = 1;
skipFrames  = 3;
numFrames   = veloReader.NumberOfFrames;
displayRate = 10;      % Update display every 100 frames
start_ins_row = 1;
firstFrame = true;
currScan = 0;
t_INS = t_INS_start;
for n = 1 : skipFrames : numFrames
    t_lidar  = PCAP_UTC_TO_GPS_TIME(seconds(veloReader.Timestamps(n))*1000000, LiDAR_CFG);
    % Check an existing time of current LiDAR
    if (t_lidar < t_INS_start)
        continue;
    elseif (t_lidar<460748)
        continue;
    elseif (t_lidar>460749)
        break;
    end
    % Time Synchronization:
    while (t_INS <= t_lidar)
        start_ins_row   = start_ins_row + 1;
        t_INS           = LLF_NED_INS(start_ins_row, 1);        
    end
    currScan = currScan+1;
    % Interpolation:
    interLocalPose = interpolatePOS(t_lidar, LLF_NED_INS, start_ins_row);
    tmpInterLocalpose                      = TF_NED2ENU(interLocalPose(:,1:10));
    LLF_ENU_INS.INT_SOL3(currScan, 1:10)    = tmpInterLocalpose;
    clearvars tmpInterLocalpose interLocalPose;
    nav_INS = getNavSolENU(LLF_ENU_INS.INT_SOL3, currScan);
    r_bn_INS = nav_INS.tran; % E,N,U
    R_n2b_INS = nav_INS.rotm;  % C_nb (n-frame: ENU to b-frame: Right-Front-Up)
    R_n2l_INS = C_bl * R_n2b_INS; % C_nl (n-frame: ENU to lidar-frame: Front-Left-Up)
    R_l2n_INS = R_n2l_INS'; % C_ln (lidar-frame: Front-Left-Up n-frame: ENU)
    r_ln_INS = (R_n2b_INS')*r_lb + r_bn_INS'; % Lidar position (from TC solution) in n-frame

    % Read point cloud
    ptCloudOrig = readFrame(veloReader,n);

    % Process point cloud
    %   - Segment and remove ground plane
    %   - Segment and remove ego vehicle
    ptCloud = helperProcessPointCloud(ptCloudOrig);

    % Visualize point cloud
    view(lidarPlayer, ptCloud);
    
    % Downsample the processed point cloud
%     ptCloud = pcdownsample(ptCloud, "random", downsamplePercent);
    
    if firstFrame
        absTform = rigidtform3d(R_l2n_INS,r_ln_INS');
        % Add first point cloud scan as a view to the view set
        % absTform = rigidtform3d(R_l2n_INS, r_ln_INS');
        vSet = addView(vSet, viewId, absTform, "PointCloud", ptCloudOrig);
        viewId = viewId + 1;
%         ptCloudPrev = pctransform(ptCloud, absTform);
        ptCloudPrev = ptCloud;
        firstFrame = false;
        continue;
    end
    
    % Use INS to estimate an initial transformation for registration
    delta_r_n = LLF_ENU_INS.INT_SOL3(currScan,2:4)'-LLF_ENU_INS.INT_SOL3(currScan-1,2:4)';
    nav_INS_prev = getNavSolENU(LLF_ENU_INS.INT_SOL3, currScan-1);
    R_n2b_INS_prev = C_bl*nav_INS_prev.rotm;  % C_nl_prev = C_bl*C_nb
    delta_r = C_bl*R_n2b_INS_prev*delta_r_n;  % under lidar frame
    R_prev_now = (C_bl*R_n2b_INS)*R_n2b_INS_prev';
    R_now_prev = R_prev_now';
    initTform = rigidtform3d(R_now_prev,delta_r);
    
    % Compute rigid transformation that registers current point cloud with
    % previous point cloud
    relTform = pcregisterndt(ptCloud, ptCloudPrev, regGridSize, ...
        "InitialTransform", initTform);
    
    % Update absolute transformation to reference frame (first point cloud)
    
    absTform = rigidtform3d( absTform.A * relTform.A );

    % Add current point cloud scan as a view to the view set
    vSet = addView(vSet, viewId, absTform, "PointCloud", ptCloudOrig);

    % Add a connection from the previous view to the current view, representing
    % the relative transformation between them
    vSet = addConnection(vSet, viewId-1, viewId, relTform);
    
    viewId = viewId + 1;
        
    ptCloudPrev = ptCloud;
    initTform   = relTform;
    
    if n>1 && mod(n, displayRate) == 1
        plot(vSet, "Parent", hAxBefore);
        drawnow update
    end
end
%% 
% The view set object |vSet|, now holds views and connections. In the Views 
% table of vSet, the |AbsolutePose| variable specifies the absolute pose of each 
% view with respect to the first view. In the |Connections| table of |vSet|, the 
% |RelativePose| variable specifies relative constraints between the connected 
% views, the |InformationMatrix| variable specifies, for each edge, the uncertainty 
% associated with a connection.

% Display the first few views and connections
head(vSet.Views)
head(vSet.Connections)
%% 
% Now, build a point cloud map using the created view set. Align the view absolute 
% poses with the point clouds in the view set using |pcalign|. Specify a grid 
% size to control the resolution of the map. The map is returned as a |pointCloud| 
% object.

ptClouds = vSet.Views.PointCloud;
absPoses = vSet.Views.AbsolutePose;
mapGridSize = 0.2;
ptCloudMap = pcalign(ptClouds, absPoses, mapGridSize);
%%
% Assuming ptCloud is your point cloud object and 'intensity' is a vector of intensity values
% minI = min(double(ptCloudMap.Intensity));
% maxI = max(double(ptCloudMap.Intensity));
% 
% % Normalize intensity to 0-1
% normI = (double(ptCloudMap.Intensity) - minI) / (maxI - minI);
% 
% % Map intensity to an RGB color (here we use a grayscale colormap, but others can be used)
% colorMap = uint8(255 * repmat(normI, 1 ,3));
% 
% % Assign color to point cloud
% ptCloudMap.Color = colorMap;
% 
% % Write to PLY file
% pcwrite(ptCloudMap, 'Odometry_vlp16_map.ply');

%% 
% Notice that the path traversed using this approach drifts over time. While 
% the path along the first loop back to the starting point seems reasonable, the 
% second loop drifts significantly from the starting point. The accumulated drift 
% results in the second loop terminating several meters away from the starting 
% point.
% 
% A map built using odometry alone is inaccurate. Display the built point cloud 
% map with the traversed path. Notice that the map and traversed path for the 
% second loop are not consistent with the first loop.

hold(hAxBefore, 'on');
pcshow(ptCloudMap);
hold(hAxBefore, 'off');

close(hAxBefore.Parent)
%% Correct Drift Using Pose Graph Optimization
% _Graph SLAM_ is a widely used technique for resolving the drift in odometry. 
% The graph SLAM approach incrementally creates a graph, where nodes correspond 
% to vehicle poses and edges represent sensor measurements constraining the connected 
% poses. Such a graph is called a _pose graph_. The pose graph contains edges 
% that encode contradictory information, due to noise or inaccuracies in measurement. 
% The nodes in the constructed graph are then optimized to find the set of vehicle 
% poses that optimally explain the measurements. This technique is called _pose 
% graph optimization_.
% 
% To create a pose graph from a view set, you can use the <docid:vision_ref#mw_46572068-cb3d-4d80-9c87-0dae0b04fb86 
% |createPoseGraph|> function. This function creates a node for each view, and 
% an edge for each connection in the view set. To optimize the pose graph, you 
% can use the <docid:nav_ref#mw_ead0d07b-dffc-4389-93ef-8fa1432c60cb |optimizePoseGraph|> 
% function.
% 
% A key aspect contributing to the effectiveness of graph SLAM in correcting 
% drift is the accurate detection of loops, that is, places that have been previously 
% visited. This is called _loop closure detection_ or _place recognition_. Adding 
% edges to the pose graph corresponding to loop closures provides a contradictory 
% measurement for the connected node poses, which can be resolved during pose 
% graph optimization.
% 
% Loop closures can be detected using descriptors that characterize the local 
% environment visible to the Lidar sensor. The _Scan Context_ descriptor [1] is 
% one such descriptor that can be computed from a point cloud using the <docid:vision_ref#mw_79434e5c-58ee-4c82-95f0-c9becef45429 
% |scanContextDescriptor|> function. This example uses a <docid:vision_ref#mw_058b7ce1-446b-4bc2-8eb1-5721a210c2dd 
% |scanContextLoopDetector|> object to manage the scan context descriptors that 
% correspond to each view. It uses the <docid:vision_ref#mw_83b4b792-2fe3-4436-b92b-3b061085a963 
% |detectLoop|> object function to detect loop closures with a two phase descriptor 
% search algorithm. In the first phase, it computes the ring key subdescriptors 
% to find potential loop candidates. In the second phase, it classifies views 
% as loop closures by thresholding the scan context distance.

% Set random seed to ensure reproducibility
rng(0);

% Create an empty view set
vSet = pcviewset;

% Create a loop closure detector
loopDetector = scanContextLoopDetector;

% Create a figure for view set display
hFigBefore = figure('Name', 'View Set Display');
hAxBefore = axes(hFigBefore);

% Initialize transformations
absTform   = rigidtform3d;  % Absolute transformation to reference frame
relTform   = rigidtform3d;  % Relative transformation between successive scans

maxTolerableRMSE  = 3; % Maximum allowed RMSE for a loop closure candidate to be accepted

viewId = 1;
skipFrames  = 3;
numFrames   = veloReader.NumberOfFrames;
displayRate = 10;      % Update display every 100 frames
start_ins_row = 1;
firstFrame = true;
currScan = 0;
t_INS = t_INS_start;
for n = 1 : skipFrames : numFrames
    t_lidar  = PCAP_UTC_TO_GPS_TIME(seconds(veloReader.Timestamps(n))*1000000, LiDAR_CFG);
    % Check an existing time of current LiDAR
    if (t_lidar < t_INS_start)
        continue;
    elseif (t_lidar<460748)
        continue;
    elseif (t_lidar>461478)
        break;
    end
    % Time Synchronization:
    while (t_INS <= t_lidar)
        start_ins_row   = start_ins_row + 1;
        t_INS           = LLF_NED_INS(start_ins_row, 1);        
    end
    currScan = currScan+1;
    % Interpolation:
    interLocalPose = interpolatePOS(t_lidar, LLF_NED_INS, start_ins_row);
    tmpInterLocalpose                      = TF_NED2ENU(interLocalPose(:,1:10));
    LLF_ENU_INS.INT_SOL3(currScan, 1:10)    = tmpInterLocalpose;
    clearvars tmpInterLocalpose interLocalPose;
    nav_INS = getNavSolENU(LLF_ENU_INS.INT_SOL3, currScan);
    r_bn_INS = nav_INS.tran; % E,N,U
    R_n2b_INS = nav_INS.rotm;  % C_nb (n-frame: ENU to b-frame: Right-Front-Up)
    R_n2l_INS = C_bl * R_n2b_INS; % C_nl (n-frame: ENU to lidar-frame: Front-Left-Up)
    R_l2n_INS = R_n2l_INS'; % C_ln (lidar-frame: Front-Left-Up n-frame: ENU)
    r_ln_INS = (R_n2b_INS')*r_lb + r_bn_INS'; % Lidar position (from TC solution) in n-frame

    % Read point cloud
    ptCloudOrig = readFrame(veloReader,n);
    
    % Process point cloud
    %   - Segment and remove ground plane
    %   - Segment and remove ego vehicle
    ptCloud = helperProcessPointCloud(ptCloudOrig);

    % Visualize point cloud
    view(lidarPlayer, ptCloud);
    
    % Downsample the processed point cloud
%     ptCloud = pcdownsample(ptCloud, "random", downsamplePercent);
    
    if firstFrame
        absTform = rigidtform3d(R_l2n_INS,r_ln_INS');
        % Add first point cloud scan as a view to the view set
        vSet = addView(vSet, viewId, absTform, "PointCloud", ptCloudOrig);
        
        % Extract the scan context descriptor from the first point cloud
        descriptor = scanContextDescriptor(ptCloudOrig);
        
        % Add the first descriptor to the loop closure detector
        addDescriptor(loopDetector, viewId, descriptor)
        
        viewId = viewId + 1;
        ptCloudPrev = ptCloud;
        firstFrame = false;
        continue;
    end
    
    % Use INS to estimate an initial transformation for registration
    delta_r_n = LLF_ENU_INS.INT_SOL3(currScan,2:4)'-LLF_ENU_INS.INT_SOL3(currScan-1,2:4)';
    nav_INS_prev = getNavSolENU(LLF_ENU_INS.INT_SOL3, currScan-1);
    R_n2b_INS_prev = C_bl*nav_INS_prev.rotm;  % C_nl_prev = C_bl*C_nb
    delta_r = C_bl*R_n2b_INS_prev*delta_r_n;  
    R_prev_now = (C_bl*R_n2b_INS)*R_n2b_INS_prev';
    R_now_prev = R_prev_now';
    initTform = rigidtform3d(R_now_prev,delta_r);
    
    % Compute rigid transformation that registers current point cloud with
    % previous point cloud
    relTform = pcregisterndt(ptCloud, ptCloudPrev, regGridSize, ...
        "InitialTransform", initTform);
    
    % Update absolute transformation to reference frame (first point cloud)
    
    absTform = rigidtform3d( absTform.A * relTform.A );
    
    % Add current point cloud scan as a view to the view set
    vSet = addView(vSet, viewId, absTform, "PointCloud", ptCloudOrig);
    
    % Add a connection from the previous view to the current view representing
    % the relative transformation between them
    vSet = addConnection(vSet, viewId-1, viewId, relTform);

    % Extract the scan context descriptor from the point cloud
    descriptor = scanContextDescriptor(ptCloudOrig);

    % Add the descriptor to the loop closure detector
    addDescriptor(loopDetector, viewId, descriptor)

    % Detect loop closure candidates
    loopViewId = detectLoop(loopDetector);
    
    % A loop candidate was found
    if ~isempty(loopViewId)
        loopViewId = loopViewId(1);
        
        % Retrieve point cloud from view set
        loopView = findView(vSet, loopViewId);
        ptCloudOrig = loopView.PointCloud;
        
        % Process point cloud
        ptCloudOld = helperProcessPointCloud(ptCloudOrig);
        
        % Downsample point cloud
%         ptCloudOld = pcdownsample(ptCloudOld, "random", downsamplePercent);
        
        % Use registration to estimate the relative pose
        [relTform, ~, rmse] = pcregisterndt(ptCloud, ptCloudOld, ...
            regGridSize, "MaxIterations", 50);
        
        acceptLoopClosure = rmse <= maxTolerableRMSE;
        if acceptLoopClosure
            % For simplicity, use a constant, small information matrix for
            % loop closure edges
            infoMat = 0.01 * eye(6);
            
            % Add a connection corresponding to a loop closure
            vSet = addConnection(vSet, loopViewId, viewId, relTform, infoMat);
        end
    end
    
    viewId = viewId + 1;
    
    ptCloudPrev = ptCloud;
    initTform   = relTform;
    
    if n>1 && mod(n, displayRate) == 1
        hG = plot(vSet, "Parent", hAxBefore);
        drawnow update
    end
end
%% 
% Create a pose graph from the view set by using the <docid:vision_ref#mw_46572068-cb3d-4d80-9c87-0dae0b04fb86 
% |createPoseGraph|> method. The pose graph is a <docid:matlab_ref#bun70tf |digraph|> 
% object with:
%
% * Nodes containing the absolute pose of each view
% * Edges containing the relative pose constraints of each connection

G = createPoseGraph(vSet);
disp(G)
% 
% In addition to the odometry connections between successive views, the view 
% set now includes loop closure connections. For example, notice the new connections 
% between the second loop traversal and the first loop traversal. These are loop 
% closure connections. These can be identified as edges in the graph whose end 
% nodes are not consecutive.

% Update axes limits to focus on loop closure connections
xlim(hAxBefore, [-50 50]);
ylim(hAxBefore, [-100 50]);

% Find and highlight loop closure connections
loopEdgeIds = find(abs(diff(G.Edges.EndNodes, 1, 2)) > 1);
highlight(hG, 'Edges', loopEdgeIds, 'EdgeColor', 'red', 'LineWidth', 3)
%% 
% Optimize the pose graph using |optimizePoseGraph|. 

optimG = optimizePoseGraph(G, 'g2o-levenberg-marquardt');

vSetOptim = updateView(vSet, optimG.Nodes);
%% 
% Display the view set with optimized poses. Notice that the detected loops 
% are now merged, resulting in a more accurate trajectory.

% plot(vSetOptim, 'Parent', hAxBefore)
plot(vSetOptim)
%% 
% The absolute poses in the optimized view set can now be used to build a more 
% accurate map. Use the <docid:vision_ref#mw_943aeeec-bf72-49f6-bffd-9a2295e6b6f4 
% |pcalign|> function to align the view set point clouds with the optimized view 
% set absolute poses into a single point cloud map. Specify a grid size to control 
% the resolution of the created point cloud map.

mapGridSize = 0.2;
ptClouds = vSetOptim.Views.PointCloud;
absPoses = vSetOptim.Views.AbsolutePose;
ptCloudMap = pcalign(ptClouds, absPoses, mapGridSize);

hFigAfter = figure('Name', 'View Set Display (after optimization)');
hAxAfter = axes(hFigAfter);
pcshow(ptCloudMap, 'Parent', hAxAfter);

% Overlay view set display
hold on
plot(vSetOptim, 'Parent', hAxAfter);

helperMakeFigurePublishFriendly(hFigAfter);
%% 
% While accuracy can still be improved, this point cloud map is significantly 
% more accurate.
%% References
%% 
% # G. Kim and A. Kim, "Scan Context: Egocentric Spatial Descriptor for Place 
% Recognition Within 3D Point Cloud Map," _2018 IEEE/RSJ International Conference 
% on Intelligent Robots and Systems (IROS)_, Madrid, 2018, pp. 4802-4809.
%% Supporting Functions and Classes
%% 
% |*helperProcessPointCloud*| processes a point cloud by removing points belonging 
% to the ground plane and the ego vehicle.

function ptCloud = helperProcessPointCloud(ptCloudIn, method)
%helperProcessPointCloud Process pointCloud to remove ground and ego vehicle
%   ptCloud = helperProcessPointCloud(ptCloudIn, method) processes 
%   ptCloudIn by removing the ground plane and the ego vehicle.
%   method can be "planefit" or "rangefloodfill".
%
%   See also pcfitplane, pointCloud/findNeighborsInRadius.

arguments
    ptCloudIn (1,1) pointCloud
    method          string      {mustBeMember(method, ["planefit","rangefloodfill"])} = "rangefloodfill"
end

isOrganized = ~ismatrix(ptCloudIn.Location);

if (method=="rangefloodfill" && isOrganized) 
    % Segment ground using floodfill on range image
    groundFixedIdx = segmentGroundFromLidarData(ptCloudIn, ...
        "ElevationAngleDelta", 11);
else
    % Segment ground as the dominant plane with reference normal
    % vector pointing in positive z-direction
    maxDistance         = 0.4;
    maxAngularDistance  = 5;
    referenceVector     = [0 0 1];

    [~, groundFixedIdx] = pcfitplane(ptCloudIn, maxDistance, ...
    referenceVector, maxAngularDistance);
end

if isOrganized
    groundFixed = false(size(ptCloudIn.Location,1),size(ptCloudIn.Location,2));
else
    groundFixed = false(ptCloudIn.Count, 1);
end
groundFixed(groundFixedIdx) = true;

% Segment ego vehicle as points within a given radius of sensor
sensorLocation = [0 0 0];
radius = 3.5;
egoFixedIdx = findNeighborsInRadius(ptCloudIn, sensorLocation, radius);

if isOrganized
    egoFixed = false(size(ptCloudIn.Location,1),size(ptCloudIn.Location,2));
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
%% 
% |helperComputeInitialEstimateFromINS*| estimates an initial transformation 
% for registration from INS readings.

%% 
% |helperMakeFigurePublishFriendly*| adjusts figures so that screenshot captured 
% by publish is correct.

function helperMakeFigurePublishFriendly(hFig)
    
if ~isempty(hFig) && isvalid(hFig)
    hFig.HandleVisibility = 'callback';
end
end
%% 
% _Copyright 2019-2021 The MathWorks, Inc._