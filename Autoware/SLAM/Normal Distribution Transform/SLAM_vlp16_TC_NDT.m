%% =============== Setup =============== %%
clc;
clear all;
close all;

%% =============== Select Input pcap File =============== %%
FileName_pcap = '2023-06-09-15-51-42_Velodyne-VLP-16-Data.pcap';
PathName_pcap = 'D:\For_thesis_student\Chin-Chia\Data\LiDAR\';
%[FileName_pcap, PathName_pcap, ~] = uigetfile('*.pcap', 'Please Select LiDAR File(.pcap)');
if isequal(FileName_pcap, 0)
    disp('User selected Cancel');
    return;
else
    disp(['User selected ', fullfile(PathName_pcap, FileName_pcap)]);
end

%% =============== Load pcap File =============== %%
[LiDAR_CFG, veloReader] = create_lidar_config(PathName_pcap, FileName_pcap);
t_lidar_start = PCAP_UTC_TO_GPS_TIME(seconds(veloReader.StartTime)*1000000, LiDAR_CFG);
t_lidar_end = PCAP_UTC_TO_GPS_TIME(seconds(veloReader.EndTime)*1000000, LiDAR_CFG);
dt = (PCAP_UTC_TO_GPS_TIME(seconds(veloReader.Timestamps(2))*1000000, LiDAR_CFG)-PCAP_UTC_TO_GPS_TIME(seconds(veloReader.Timestamps(1))*1000000, LiDAR_CFG));
lidarRate = round(1/dt);

if isequal(veloReader.DeviceModel, 'HDL64E')
    % ===== Define Vertical and Horizontal Resolution ===== %
    LiDAR_CFG.vResolution = 64;
    LiDAR_CFG.hResolution = 1024;

    % ===== Define Vertical and Horizontal Field-of-View ===== %
    LiDAR_CFG.vFoVUp = 2;     
    LiDAR_CFG.vFoVDown = -24.9; 
    LiDAR_CFG.vFoV = [LiDAR_CFG.vFoVUp LiDAR_CFG.vFoVDown];
    LiDAR_CFG.hFoV = 360;

elseif isequal(veloReader.DeviceModel, 'VLP16')
    % ===== Define Vertical and Horizontal Resolution ===== %
    LiDAR_CFG.vResolution = 16;
    LiDAR_CFG.hResolution = 1024;

    % ===== Define Vertical and Horizontal Field-of-View ===== %
    LiDAR_CFG.vFoVUp = 15;     
    LiDAR_CFG.vFoVDown = -15; 
    LiDAR_CFG.vFoV = [LiDAR_CFG.vFoVUp LiDAR_CFG.vFoVDown];
    LiDAR_CFG.hFoV = 360;

end

%% =============== Select Input Integrated Navigation File =============== %%
FileName_IE = 'REF_1.txt';
PathName_IE = 'D:\For_thesis_student\Chin-Chia\Data\POS[IE]\';
%[FileName_IE, PathName_IE, ~] = uigetfile('*.txt', 'Please Select Reference File(.txt)');
if isequal(FileName_IE, 0)
    disp('User selected Cancel');
    return;
else
    disp(['User selected ', fullfile(PathName_IE, FileName_IE)]);
end

%% =============== Load Integrated Navigation Data =============== %%
initPOS = [23.0028474266, 120.214659453];
initPOS4WGS84 = [initPOS(1), initPOS(2), 38.070];
initPOS4TWD97 = wgs84_to_twd97(initPOS);

File_INS = load([PathName_IE, FileName_IE]);
t_INS_start = File_INS(1, 1);
LLF_NED_INS = wgs2llf(File_INS, initPOS4WGS84);

for i = 1:size(LLF_NED_INS, 1)
    if (LLF_NED_INS(i, 10) < 0)
        LLF_NED_INS(i, 10) = LLF_NED_INS(i, 10) + 360;
    end
    LLF_NED_INS(i, 2) = LLF_NED_INS(i, 2) + initPOS4TWD97(2);
    LLF_NED_INS(i, 3) = LLF_NED_INS(i, 3) + initPOS4TWD97(1);
end

%% =============== Initial LiDAR Mounting Parameters =============== %%
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

%% =============== Select Input pcd File =============== %%
%FileName_initpcd = 'vlp16_DG_14850_15170.pcd';
%PathName_initpcd = 'D:\For_thesis_student\Chin-Chia\Normal Distribution Transform\input-pcd\';
[FileName_initpcd, PathName_initpcd, ~] = uigetfile('*.pcd', 'Please Select Initial Map File(.pcd)');
if isequal(FileName_initpcd, 0)
    disp('User selected Cancel');
    return;
else
    disp(['User selected ', fullfile(PathName_initpcd, FileName_initpcd)]);
end

%% =============== Load pcd File =============== %%
ptCloud_initpcd = pcread(fullfile(PathName_initpcd, FileName_initpcd));
%figure;pcshow(ptCloud_initpcd);

%% =============== Select Output pcd Folder =============== %%
PathName_output = 'D:\For_thesis_student\Chin-Chia\Normal Distribution Transform\output-pcd';
%PathName_output = uigetdir(addpath(genpath(pwd)), 'Please Select Output Folder');
if isequal(PathName_output, 0)
    disp('User selected Cancel');
    return;
else
    disp(['User selected ', PathName_output]);
end

%% =============== Mapping Process =============== %%
xlimits = [-15, 15];
ylimits = [-15, 15];
zlimits = [-10, 20];
lidarPlayer = pcplayer(xlimits, ylimits, zlimits);

xlabel(lidarPlayer.Axes, 'X (m)');
ylabel(lidarPlayer.Axes, 'Y (m)');
zlabel(lidarPlayer.Axes, 'Z (m)');

% ===== Create an Empty View Set ===== %
vSet = pcviewset;

% ===== Create a Figure for View Set Display ===== %
hFigBefore = figure('Name', 'View Set Display');
hAxBefore = axes(hFigBefore);

% ===== Initialize Transformations ===== %
absTform = rigidtform3d;  % Absolute transformation to reference frame
relTform = rigidtform3d;  % Relative transformation between successive scans

viewId = 1;
stratFrame = 09740;
endFrame = 18910;
skipFrame = 1;
numFrame = veloReader.NumberOfFrames;
displayRate = 100;
firstFrame = true;
start_ins_row = 1;
t_INS = t_INS_start;
currScan = 0;

for n = stratFrame: skipFrame: endFrame
    % ===== Select Process Period ===== %
    if ~((n >= 9740 && n <= 12560) || (n >= 14020 && n <= 14740) || (n >= 15830 && n <= 17210) || (n >= 17380 && n <= 17900) || (n >= 18030 && n <= 18910))
        continue;
    end
    
%     if ~((n >= 15830 && n <= 17210))
%         continue;
%     end

    % ===== Align LiDAR Time and INS Time ===== %
    t_lidar = PCAP_UTC_TO_GPS_TIME(seconds(veloReader.Timestamps(n))*1000000, LiDAR_CFG);
    title(lidarPlayer.Axes, sprintf('Lidar Sensor Data at Frame NO.%5d', n));
    
    while (t_INS <= t_lidar)
        start_ins_row = start_ins_row + 1;
        t_INS = LLF_NED_INS(start_ins_row, 1);
    end
    currScan = currScan + 1;

    % ===== Interpolation ===== %
    interLocalPose = interpolatePOS(t_lidar, LLF_NED_INS, start_ins_row);
    tmpInterLocalpose = TF_NED2ENU(interLocalPose(:, 1: 10));
    LLF_ENU_INS.INT_SOL3(currScan, 1: 10) = tmpInterLocalpose;
    clearvars tmpInterLocalpose interLocalPose;
    nav_INS = getNavSolENU(LLF_ENU_INS.INT_SOL3, currScan);
    r_bn_INS = nav_INS.tran; % E,N,U
    R_n2b_INS = nav_INS.rotm;  % C_nb (n-frame: ENU to b-frame: Right-Front-Up)
    R_n2l_INS = C_bl * R_n2b_INS; % C_nl (n-frame: ENU to lidar-frame: Front-Left-Up)
    R_l2n_INS = R_n2l_INS'; % C_ln (lidar-frame: Front-Left-Up n-frame: ENU)
    r_ln_INS = (R_n2b_INS')*r_lb + r_bn_INS'; % Lidar position (from TC solution) in n-frame
    
    % ===== Initial Point Cloud ===== %
    ptCloudObj = readFrame(veloReader, n);

    % ===== Process Point Cloud ===== %
    ptCloud = helperProcessPointCloud_SetRange(ptCloudObj);

    % ===== Visualize the Point Cloud ===== %
    view(lidarPlayer, ptCloud.Location, ptCloud.Intensity);
    
    % ===== Setup Initial Value ===== %
    if firstFrame
        absTform = rigidtform3d(R_l2n_INS, r_ln_INS');
        vSet = addView(vSet, viewId, absTform, "PointCloud", ptCloudObj);
        viewId = viewId + 1;
        ptCloudPrev = ptCloud;
        firstFrame = false;
        continue;
    end

    % ===== Use INS to Estimate an Initial Transformation for Registration ===== %
    delta_r_n = LLF_ENU_INS.INT_SOL3(currScan, 2: 4)' - LLF_ENU_INS.INT_SOL3(currScan-1, 2: 4)';
    nav_INS_prev = getNavSolENU(LLF_ENU_INS.INT_SOL3, currScan-1);
    R_n2l_INS_prev = C_bl * nav_INS_prev.rotm;  % C_nl_prev = C_bl*C_nb
    delta_r = R_n2l_INS_prev * delta_r_n;  % under lidar frame
    R_prev_now = (C_bl*R_n2b_INS) * R_n2l_INS_prev';
    R_now_prev = R_prev_now';
    initTform = rigidtform3d(R_now_prev, delta_r);

    % ===== Compute Rigid Transformation that Registers Current Point Cloud with Previous Point Cloud ===== %
    relTform = pcregisterndt(ptCloud, ptCloud_initpcd, 25);

    % ===== Update Absolute Transformation to Reference Frame (first point cloud) ===== %
    absTform = rigidtform3d(R_l2n_INS, r_ln_INS');

    % ===== Add Current Point Cloud Scan as a View to the View Set ===== %
    vSet = addView(vSet, viewId, absTform, "PointCloud", ptCloudObj);

    % ===== Add a Connection from the Previous View to the Current View, R epresenting the Relative Transformation Between Them ===== %
    vSet = addConnection(vSet, viewId-1, viewId, relTform);
    viewId = viewId + 1;
        
    ptCloudPrev = ptCloud;
    initTform = relTform;

    % ===== Update Rate ===== %
    if n > 1 && mod(n, displayRate) == 0
        plot(vSet, "Parent", hAxBefore);
        drawnow update
    end

    pause(0.01);
end

%% =============== Save Output pcd =============== %%
ptClouds = vSet.Views.PointCloud;
absPoses = vSet.Views.AbsolutePose;
mapGridSize = 0.2;
ptCloudMap = pcalign(ptClouds, absPoses, mapGridSize);

FileName_output = sprintf('vlp16_NDT_%05d_%05d_.pcd', stratFrame, endFrame);
file_output = [PathName_output, '\', FileName_output];
output_pcd(ptCloudMap, file_output);
