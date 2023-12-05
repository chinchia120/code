%============================================================================================================================
% Purposes: Multi-LiDAR Calibration Using HD Point Cloud Map-based Least-Squares Adjustment
% Developed & Coded by: Surachet Srinara
% Date: June 4, 2022 - 02:33AM
% Last modified: June 12, 2022 - 05:45PM
% Copyright@Chet2022
%============================================================================================================================ 
fclose('all'); clearvars; close all; clc
addpath(genpath(pwd))

%% Load LiDAR Data:
% == First LiDAR: VLP-16 (Front-LiDAR)
lidarFile1 = '2023-06-09-15-51-42_Velodyne-VLP-16-Data.pcap';
lidarPath1 = 'C:\Users\user\Documents\NCKU-Data\competition\SG\data\';
%[lidarFile1, lidarPath1, ~]     = uigetfile('*.pcap', 'Please select your first LiDAR file');
if lidarFile1 == 0
    return;
end
[LiDAR_CFG1, veloReader1]       = create_lidar_config(lidarPath1, lidarFile1);
t_lidar_start1                  = PCAP_UTC_TO_GPS_TIME(seconds(veloReader1.StartTime)*1000000, LiDAR_CFG1); 
t_lidar_end1                    = PCAP_UTC_TO_GPS_TIME(seconds(veloReader1.EndTime)*1000000, LiDAR_CFG1);
dt1                             = (PCAP_UTC_TO_GPS_TIME(seconds(veloReader1.Timestamps(2))*1000000, LiDAR_CFG1)-...
                                   PCAP_UTC_TO_GPS_TIME(seconds(veloReader1.Timestamps(1))*1000000, LiDAR_CFG1));
lidarRate1                      = round(1/dt1);
if isequal(veloReader1.DeviceModel, 'HDL64E')
    % Define vertical and horizontal resolution.
    LiDAR_CFG1.vResolution = 64;
    LiDAR_CFG1.hResolution = 1024;
    % Define vertical and horizontal field-of-view.
    LiDAR_CFG1.vFoVUp      = 2;     
    LiDAR_CFG1.vFoVDown    = -24.9; 
    LiDAR_CFG1.vFoV        = [LiDAR_CFG1.vFoVUp LiDAR_CFG1.vFoVDown];
    LiDAR_CFG1.hFoV        = 360;
elseif isequal(veloReader1.DeviceModel, 'VLP16')
    % Define vertical and horizontal resolution.
    LiDAR_CFG1.vResolution = 16;
    LiDAR_CFG1.hResolution = 1024;
    % Define vertical and horizontal field-of-view.
    LiDAR_CFG1.vFoVUp      = 15;     
    LiDAR_CFG1.vFoVDown    = -15; 
    LiDAR_CFG1.vFoV        = [LiDAR_CFG1.vFoVUp LiDAR_CFG1.vFoVDown];
    LiDAR_CFG1.hFoV        = 360;
end

%% Calculate LiDAR Time Span:
% t_lidars                        = [t_lidar_start1, t_lidar_end1; t_lidar_start2, t_lidar_end2; t_lidar_start3, t_lidar_end3];
t_lidars                        = [t_lidar_start1, t_lidar_end1];
t_lidars_start                  = max(t_lidars(:,1));
t_lidars_end                    = min(t_lidars(:,2));

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
% la10    = [0.163; 1.315; 0.784];
% la20    = [0.163; 0.545; 0.784];
la    = [0.163; 1.095; (1.231+1.158)/2];
% r1_lb   = la10;
% r2_lb   = la20;
r_lb   = la;

% == Bore-sight angles: [bs_right; bs_fwd; bs_up]
bs    = [0; 0; 0];
bsRotm = eul2rotmENU(bs);
C_bl_intial_guess  = bsRotm.matrix;

%% LiDAR Configuration:
% == Specify time slots:                                                       
startTime       = 462416;                                                                                              
endTime         = 462426;                                                           
freqFrame       = 0.1;
timeSlots       = startTime:freqFrame:endTime;                                   

% Specify scoping ranges: 
minRange        = 3.5;                                                          
maxRange       = 70;                                                               

% Specify voxel size for NDT:
LiDAR_CFG1.VoxelSizeNDT = 1;


%% Load HD Point Cloud Map:
%fnameHD = 'Odometry_vlp16_map.pcd'
%pnameHD = 'C:\Users\user\Documents\NCKU-Data\competition\SG\data\';
[fnameHD, pnameHD] = uigetfile({'*.pcd'}, 'Please select your HD map file','MultiSelect','on');
% for i=1:size(fnameHD,2)
%     ptCloudOut(i,1) = pcread(fullfile(pnameHD,fnameHD{i}));
% end
% ptCloudOut = pccat(ptCloudOut);
% pcwrite(ptCloudOut,fullfile(pnameHD,"HDmap_transform.pcd"))
ptCloudHD = pcread(fullfile(pnameHD,fnameHD));
clearvars pnameHD fnameHD;

%%  Load Integrated Navigation Data: TC-INS/GNSS (IE)
initPOS                         = [23.0028474266, 120.214659453];
initPOS4WGS84                   = [initPOS(1), initPOS(2), 38.070];
[FileNameGT, PathNameGT, ~]     = uigetfile('*.txt', 'Please select your reference file (IE)');
f_ref                           = load([PathNameGT FileNameGT]);
t_gt                            = f_ref(1, 1);
LLF_NED_GT0                     = wgs2llf(f_ref, initPOS4WGS84);
LLF_NED_GT                      = LLF_NED_GT0;
for i = 1:size(LLF_NED_GT,1)
    if (LLF_NED_GT(i,10) < 0)
        LLF_NED_GT(i,10) = LLF_NED_GT(i,10) + 360;
    else
        continue;
    end
end
clearvars LLF_NED_GT0 i;

%% Output Configuration:
formatOut       = 'yymmdd';
headerFilename  = datestr(now,formatOut);
% insDir          = strcat(PathNameGT, headerFilename, '_INS_NAV', '.txt');
% ndtDir          = strcat(PathNameGT, headerFilename, '_NDT_NAV', '.txt');
% fid_ins         = fopen(insDir, 'w');
% fid_ndt         = fopen(ndtDir, 'w');
% outStr          = '%.8f\t%.20f\t%.20f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\n';

%% Main Process: 
currScan       = 0; 
cntMatch      = 0;
t_gt1           = t_gt;
start_gt_row   = 1;
% cntTimeSlot     = 1;

while (hasFrame(veloReader1)) % VLP16
    % Check an existing scan frame
    if ((currScan+1) > veloReader1.NumberOfFrames)
        fprintf('End of Frame!');
        break;
    end
    currScan   = currScan + 1;
    t_lidar    = PCAP_UTC_TO_GPS_TIME(seconds(veloReader1.Timestamps(currScan))*1000000, LiDAR_CFG1);
    fprintf('Primary LiDAR Frame No.%0.0f [%0.6f]\n', currScan, t_lidar);
    % Check an existing time of current LiDAR
    if (t_lidar < f_ref(1,1))
        continue;
    end	
	% Time Synchronization:
    while (t_gt1 <= t_lidar)
        start_gt_row   = start_gt_row + 1;
        t_gt1           = LLF_NED_GT(start_gt_row, 1);        
    end
    % Interpolation:
%     interWGSPose                           = interpolatePOS(t_lidar, f_ref, start_gt_row);
    interLocalPose                         = interpolatePOS(t_lidar, LLF_NED_GT, start_gt_row);
    tmpInterLocalpose                      = TF_NED2ENU(interLocalPose(:,1:10));
    LLF_ENU_GT.INT_SOL3(currScan, 1:10)    = tmpInterLocalpose;
    clearvars tmpInterLocalpose interLocalPose;
    % Get an initial navigation data:     
    nav_gt      = getNavSolENU(LLF_ENU_GT.INT_SOL3, currScan);
    r_bn_gt     = nav_gt.tran; % E,N,U
    v_bn_gt     = nav_gt.vel;  % VE, VN, VU
%     a_bn_gt     = nav_gt.att'; % roll, pitch, heading
    R_n2b_gt    = nav_gt.rotm;  % C_nb (n-frame: ENU to b-frame: Right-Front-Up)
    R_n2l_gt    = C_bl_intial_guess * R_n2b_gt; % C_nl (n-frame: ENU to lidar-frame: Front-Left-Up)
    R_l2n_gt    = R_n2l_gt'; % C_ln (lidar-frame: Front-Left-Up n-frame: ENU)
    a_ln_gt     = Cbn2eul(R_l2n_gt); 
    r_ln_gt     = (R_n2b_gt')*r_lb + r_bn_gt'; % Lidar position (from TC solution) in n-frame

    % If it is in the time slot of calibration
    if any(((t_lidar+0.01) >= startTime) & (t_lidar <= endTime))
        cntMatch = cntMatch + 1;
        fprintf('-----------------LiDARs-Time Synchronization-----------------\n');
        fprintf('                              VLP16-front\n');    
        fprintf('    Scan No.          %0.0f\n', currScan);
        fprintf('    Timestamp      %0.3f\n', t_lidar);

        % Load Raw Point Cloud Object:
        ptCloudRaw  = readFrame(veloReader1,currScan);

        % Organize Raw Point Clouds:
        ptCloudOrg  = ptCloudOrgConversion(ptCloudRaw, LiDAR_CFG1);

        % Direct Geo-reference Point Cloud:
        ptCloudDG          = pcalign(ptCloudOrg, rigidtform3d(R_l2n_gt, r_ln_gt'));
        ROI3                = [ptCloudDG.XLimits(1), ptCloudDG.XLimits(2), ...
                                               ptCloudDG.YLimits(1), ptCloudDG.YLimits(2), ...
                                               ptCloudHD.ZLimits(1), ptCloudHD.ZLimits(2)];
        ind3                = findPointsInROI(ptCloudHD, ROI3);
        extractedHDmap      = select(ptCloudHD, ind3);
        clearvars  ind3 ROI3

        % Preprocess LiDAR
        ROI3                = [extractedHDmap.XLimits(1), extractedHDmap.XLimits(2), ...
                                               extractedHDmap.YLimits(1), extractedHDmap.YLimits(2), ...
                                               extractedHDmap.ZLimits(1), extractedHDmap.ZLimits(2)];
        ind3                = findPointsInROI(ptCloudDG, ROI3);
        ptCloudDG      = select(ptCloudDG, ind3);
        ptCloud          = pctransform(ptCloudDG, invert(rigidtform3d(R_l2n_gt, r_ln_gt')));
        clearvars  ind3 ROI3

        % LiDAR Matching
        initTform  = rigidtform3d(R_l2n_gt, r_ln_gt');
        % Process point cloud
        %   - Segment and remove ground plane
        %   - Segment and remove ego vehicle
        ptCloud = helperProcessPointCloud(ptCloud);

        [tform,ndt_pc,rmse] = pcregisterndt(ptCloud,extractedHDmap,LiDAR_CFG1.VoxelSizeNDT,"InitialTransform",initTform,'MaxIterations',50, 'Tolerance',  [0.01, 0.1]);
        r_ln_ndt = tform.Translation';
        R_l2n_ndt = tform.R; % After NDT matching: Cln
        R_n2l_ndt = R_l2n_ndt'; % Cnl
        a_ln_ndt    = Cbn2eul(R_l2n_ndt);
        R_b2l_ndt = R_n2l_ndt*(R_n2b_gt'); % C_bl = C_nl * C_bn
        a_lb_ndt = getEulerFromRotm(R_b2l_ndt, [0;0;0]);

        r_lb_ndt(:,cntMatch) = R_n2b_gt*(r_ln_ndt-r_bn_gt');
        a_lb(:,cntMatch) = a_lb_ndt;
        match_quality(cntMatch) = rmse;
        r_ln_leverarm     = (R_n2b_gt')*r_lb_ndt(:,cntMatch) + r_bn_gt';
        ptCloudDG_ndt  = pcalign(ptCloud, rigidtform3d(R_l2n_ndt, r_ln_leverarm'));
        fprintf('\n-------------------- Calibration result --------------------\n');
        fprintf('\nLever-arms (m) : %0.5f (R), %0.5f (F), %0.5f (U)\n', r_lb_ndt(:,cntMatch)');
        fprintf('\nBore-sight angles (deg.):  [%0.3f; %0.3f; %0.3f]\n', a_lb(:,cntMatch)');
        fprintf('\nrmse (m):  %0.3f\n', rmse);
        fprintf('--------------------------------------------------------------\n');
    else
       continue;
    end
end
save vlp16.mat r_lb_ndt a_lb;
%% average the result
r_lb_calib(1,1) = mean(r_lb_ndt(1,r_lb_ndt(1,:)>0));
r_lb_calib(2,1) = mean(r_lb_ndt(2,:));
r_lb_calib(3,1) = mean(r_lb_ndt(3,:));
a_lb_calib = mean(a_lb,2);