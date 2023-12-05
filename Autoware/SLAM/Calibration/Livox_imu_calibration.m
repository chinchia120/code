%============================================================================================================================
% Purposes: Multi-LiDAR Calibration Using HD Point Cloud Map-based Least-Squares Adjustment
% Developed & Coded by: Surachet Srinara 
% Modified by : Daniel Lu
% Date: May 31, 2023
% Copyright@Chet2022
%============================================================================================================================ 
fclose('all'); clearvars; close all; clc
addpath(genpath(pwd))
%% Load LiDAR Data:
selpath = uigetdir(pwd,'Get livox pointcloud');
if isequal(selpath,0)
    disp('User selected Cancel');
else
    disp(['User selected ', selpath]);
end
livox_pcds = datastore(selpath,"Type","file","ReadFcn",@load_livox,"IncludeSubfolders",true,"FileExtensions",'.pcd');
num_frame = size(livox_pcds.Files,1);
reset(livox_pcds);
t_end=string(livox_pcds.Files{end});
t_end = split(t_end,'.');
t_end = split(t_end(1),'\');
t_end = str2double(t_end(end))*10^-9;
utc = unix2time(t_end);
gpst = utc2gpst(utc);
[~, t_end] = time2gpst(gpst);
t_start=string(livox_pcds.Files{1});
t_start = split(t_start,'.');
t_start = split(t_start(1),'\');
t_start = str2double(t_start(end))*10^-9;
utc = unix2time(t_start);
gpst = utc2gpst(utc);
[week, t_start] = time2gpst(gpst);
clear gpst utc;
%% Calculate LiDAR Time Span:
t_lidars                        = [t_start, t_end];
t_lidars_start                  = max(t_lidars(:,1));
t_lidars_end                    = min(t_lidars(:,2));
disp(['Week: ', num2str(week), ' Elapsed time: ',num2str(t_end-t_start), ' seconds']);
clearvars t_end t_start week;
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
la    = [0.18;1.441;0.967];
% la    = [0.12;1.2561;0.9628];
r_lb   = la;

% == Bore-sight angles: [bs_right; bs_fwd; bs_up]
bs    = [0; 0; 90];
% bs = [-0.082,3.1703,88.4];
bsRotm = eul2rotmENU(bs);
C_bl  = bsRotm.matrix;
%% LiDAR Configuration:
% == Specify time slots:                                                       
startTime       = 490502;                                                                                              
endTime         = 490936;                                                           
freqFrame       = 0.1;
timeSlots       = startTime:freqFrame:endTime; 
% Specify scoping ranges: 
minRange        = 3.5;                                                          
maxRange1       = 70;                                                          
% Specify voxel size for NDT:
LiDAR_CFG.VoxelSizeNDT = 1;
clearvars freqFrame;
%%  Load Integrated Navigation Data: TC-INS/GNSS (IE)
initPOS                         = [23.0028474266, 120.214659453];
initPOS4WGS84                   = [initPOS(1), initPOS(2), 38.070];
[FileNameGT, PathNameGT, ~]     = uigetfile('*.txt', 'Please select your reference file (IE)');
f_ref                           = load([PathNameGT FileNameGT]);
t_gt                            = f_ref(1, 1);
% initPOS4WGS84 = f_ref(1,2:4);
LLF_NED_GT0                     = wgs2llf(f_ref, initPOS4WGS84);
LLF_NED_GT                      = LLF_NED_GT0;
for i = 1:size(LLF_NED_GT,1)
    if (LLF_NED_GT(i,10)) < 0
        LLF_NED_GT(i,10) = LLF_NED_GT(i,10) + 360;
    else
        continue;
    end
end
clearvars LLF_NED_GT0 i;
%% Load HD Point Cloud Map:
[fnameHD, pnameHD] = uigetfile({'*.pcd'}, 'Please select your HD map file','MultiSelect','on');
% for i=1:size(fnameHD,2)
%     ptCloudOut(i,1) = pcread(fullfile(pnameHD,fnameHD{i}));
% end
% ptCloudOut = pccat(ptCloudOut);
% pcwrite(ptCloudOut,fullfile(pnameHD,"HDmap_transform.pcd"))
ptCloudHD = pcread(fullfile(pnameHD,fnameHD));
clearvars pnameHD fnameHD;
%% Main Process: 
currScan       = 0; 
cntMatch      = 0;
t_gt1           = t_gt;
start_gt_row   = 1;
cntTimeSlot     = 1;

reset(livox_pcds);
while(hasdata(livox_pcds))
    % Check an existing scan frame
    if (currScan > num_frame)
        disp('End of Frame!');
        break;
    end
    ptCloud = read(livox_pcds);
    t_lidar=string(livox_pcds.Files{currScan});
    t_lidar = split(t_lidar,'.');
    t_lidar = split(t_lidar(1),'\');
    t_lidar = str2double(t_lidar(end))*10^-9;
    utc = unix2time(t_lidar);
    gpst = utc2gpst(utc);
    [t_week, t_lidar] = time2gpst(gpst);
    fprintf('Primary LiDAR Frame No.%0.0f [%0.6f]\n', currScan, t_lidar);
    currScan   = currScan + 1;
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
    interWGSPose                          = interpolatePOS(t_lidar, f_ref, start_gt_row);
    interLocalPose                         = interpolatePOS(t_lidar, LLF_NED_GT, start_gt_row);
    tmpInterLocalpose                      = TF_NED2ENU(interLocalPose(:,1:10));
    LLF_ENU_GT.INT_SOL3(currScan, 1:10)    = tmpInterLocalpose;
    clearvars tmpInterLocalpose interLocalPose;
    % Get an initial navigation data:     
    nav_gt      = getNavSolENU(LLF_ENU_GT.INT_SOL3, currScan);
    r_bn_gt     = nav_gt.tran; % E,N,U
    v_bn_gt     = nav_gt.vel;  % VE, VN, VU
    a_bn_gt     = nav_gt.att'; % roll, pitch, -heading
    R_n2b_gt    = nav_gt.rotm;  % C_nb (n-frame: ENU to b-frame: Right-Front-Up)
    R_n2l_gt    = C_bl * R_n2b_gt; % C_nl (n-frame: ENU to lidar-frame: Front-Left-Up)
    R_l2n_gt    = R_n2l_gt'; % C_ln (lidar-frame: Front-Left-Up n-frame: ENU)
    a_ln_gt     = Cbn2eul(R_l2n_gt); 
    r_ln_gt     = (R_n2b_gt')*r_lb + r_bn_gt'; % Lidar position (from TC solution) in n-frame

    % If it is in the time slot of calibration
    if any(((t_lidar+0.01) >= startTime) & (t_lidar <= endTime))
        cntMatch = cntMatch + 1;
        fprintf('-----------------LiDARs-Time Synchronization-----------------\n');
        fprintf('                              Livox-front\n');    
        fprintf('    Scan No.          %0.0f\n', currScan);
        fprintf('    Timestamp      %0.3f\n', t_lidar);

        % Direct Geo-reference Point Cloud:
        ptCloudDG          = pcalign(ptCloud, rigidtform3d(R_l2n_gt, r_ln_gt'));
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
%         initTform  = rigidtform3d(eye(3), [0,0,0]);
        % Process point cloud
        %   - Segment and remove ground plane
        %   - Segment and remove ego vehicle
        ptCloud = helperProcessPointCloud(ptCloud);
        [tform,ndt_pc,rmse] = pcregisterndt(ptCloud,extractedHDmap,LiDAR_CFG.VoxelSizeNDT,"InitialTransform",initTform,'MaxIterations',50);
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
save livox.mat r_lb_ndt a_lb;
%% average the result
r_lb_calib(1,1) = mean(r_lb_ndt(1,r_lb_ndt(1,:)>0));
r_lb_calib(2,1) = mean(r_lb_ndt(2,:));
r_lb_calib(3,1) = mean(r_lb_ndt(3,:));
a_lb_calib = mean(a_lb,2);
%% DG based on calibration result
a_lb_calib =  [-0.082,3.1703,88.4];
r_lb_calib = [0.12;1.2561;0.9628];
bsRotm = eul2rotmENU(a_lb_calib);
C_bl_calib  = bsRotm.matrix;
R_n2l_calib    = C_bl_calib * R_n2b_gt; % C_nl (n-frame: ENU to lidar-frame: Front-Left-Up)
R_l2n_calib    = R_n2l_calib'; % C_ln (lidar-frame: Front-Left-Up n-frame: ENU)
r_ln_calib     = (R_n2b_gt')*r_lb_calib + r_bn_gt'; % Lidar position (from TC solution) in n-frame
ptCloudDG_calib  = pcalign(ptCloud, rigidtform3d(R_l2n_calib, r_ln_calib'));
