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
[lidarFile1, lidarPath1, ~]     = uigetfile('*.pcap', 'Please select your first LiDAR file');
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

% == Second LiDAR: VLP-16 (Back-LiDAR) 
% [lidarFile2, lidarPath2, ~]     = uigetfile('*.pcap', 'Please select your second LiDAR file');
% if lidarFile2 == 0
%     return;
% end
% [LiDAR_CFG2, veloReader2]       = create_lidar_config(lidarPath2, lidarFile2);
% t_lidar_start2                  = PCAP_UTC_TO_GPS_TIME(seconds(veloReader2.StartTime)*1000000, LiDAR_CFG2); 
% t_lidar_end2                    = PCAP_UTC_TO_GPS_TIME(seconds(veloReader2.EndTime)*1000000, LiDAR_CFG2);
% dt2                             = (PCAP_UTC_TO_GPS_TIME(seconds(veloReader2.Timestamps(2))*1000000, LiDAR_CFG2)-...
%                                    PCAP_UTC_TO_GPS_TIME(seconds(veloReader2.Timestamps(1))*1000000, LiDAR_CFG2));
% lidarRate2                      = round(1/dt2);
% if isequal(veloReader2.DeviceModel, 'HDL64E')
%     % Define vertical and horizontal resolution.
%     LiDAR_CFG2.vResolution = 64;
%     LiDAR_CFG2.hResolution = 1024;
%     % Define vertical and horizontal field-of-view.
%     LiDAR_CFG2.vFoVUp      = 2;     
%     LiDAR_CFG2.vFoVDown    = -24.9; 
%     LiDAR_CFG2.vFoV        = [LiDAR_CFG2.vFoVUp LiDAR_CFG2.vFoVDown];
%     LiDAR_CFG2.hFoV        = 360;
% elseif isequal(veloReader2.DeviceModel, 'VLP16')
%     % Define vertical and horizontal resolution.
%     LiDAR_CFG2.vResolution = 16;
%     LiDAR_CFG2.hResolution = 1024;
%     % Define vertical and horizontal field-of-view.
%     LiDAR_CFG2.vFoVUp      = 15;     
%     LiDAR_CFG2.vFoVDown    = -15; 
%     LiDAR_CFG2.vFoV        = [LiDAR_CFG2.vFoVUp LiDAR_CFG2.vFoVDown];
%     LiDAR_CFG2.hFoV        = 360;
% end
% 
% == Third LiDAR: HDL-64E (Middle-LiDAR)
% [lidarFile3, lidarPath3, ~]     = uigetfile('*.pcap', 'Please select your third LiDAR file');
% if lidarFile3 == 0
%     return;
% end
% [LiDAR_CFG3, veloReader3]       = create_lidar_config(lidarPath3, lidarFile3);
% t_lidar_start3                  = PCAP_UTC_TO_GPS_TIME(seconds(veloReader3.StartTime)*1000000, LiDAR_CFG3); 
% t_lidar_end3                    = PCAP_UTC_TO_GPS_TIME(seconds(veloReader3.EndTime)*1000000, LiDAR_CFG3);
% dt3                             = (PCAP_UTC_TO_GPS_TIME(seconds(veloReader3.Timestamps(2))*1000000, LiDAR_CFG3)-...
%                                    PCAP_UTC_TO_GPS_TIME(seconds(veloReader3.Timestamps(1))*1000000, LiDAR_CFG3));
% lidarRate3                      = round(1/dt3);
% if isequal(veloReader3.DeviceModel, 'HDL64E')
%     % Define vertical and horizontal resolution.
%     LiDAR_CFG3.vResolution = 64;
%     LiDAR_CFG3.hResolution = 1024;
%     % Define vertical and horizontal field-of-view.
%     LiDAR_CFG3.vFoVUp      = 2;     
%     LiDAR_CFG3.vFoVDown    = -24.9; 
%     LiDAR_CFG3.vFoV        = [LiDAR_CFG3.vFoVUp LiDAR_CFG3.vFoVDown];
%     LiDAR_CFG3.hFoV        = 360;
% elseif isequal(veloReader3.DeviceModel, 'VLP16')
%     % Define vertical and horizontal resolution.
%     LiDAR_CFG3.vResolution = 16;
%     LiDAR_CFG3.hResolution = 1024;
%     % Define vertical and horizontal field-of-view.
%     LiDAR_CFG3.vFoVUp      = 15;     
%     LiDAR_CFG3.vFoVDown    = -15; 
%     LiDAR_CFG3.vFoV        = [LiDAR_CFG3.vFoVUp LiDAR_CFG3.vFoVDown];
%     LiDAR_CFG3.hFoV        = 360;
% end

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
la30    = [0.163; 1.095; (1.231+1.158)/2];
% r1_lb   = la10;
% r2_lb   = la20;
r3_lb   = la30;

% == Bore-sight angles: [bs_right; bs_fwd; bs_up]
% bs10    = [0; 0; 0];
% bs20    = [0; 0; -180];
bs30    = [0; 0; 0];
% bsRotm1 = eul2rotmENU(bs10);
% bsRotm2 = eul2rotmENU(bs20);
bsRotm3 = eul2rotmENU(bs30);
% R1_b2l  = bsRotm1.matrix;
% R2_b2l  = bsRotm2.matrix;
R3_b2l  = bsRotm3.matrix;

%% LiDAR Configuration:
% == Specify time slots:                                                       
startTime       = 102791;                                                                                              
endTime         = 102978;                                                           
freqFrame       = 0.1;
timeSlots       = startTime:freqFrame:endTime;                                   

% Specify scoping ranges: 
minRange        = 3.5;                                                          
maxRange1       = 70;                                                           
% maxRange2       = 70;                                                           
% maxRange3       = 90;    

% Specify voxel size for NDT:
LiDAR_CFG1.VoxelSizeNDT = 1;
% LiDAR_CFG2.VoxelSizeNDT = 1;
% LiDAR_CFG3.VoxelSizeNDT = 1;

%% Load HD Point Cloud Map:
[fnameHD, pnameHD] = uigetfile({'*.mat'}, 'Please select your HD map file');
if fnameHD == 0
    return;
end
load([pnameHD fnameHD]);
% xlimits     = [adjustedHDmap.XLimits(1)-50 adjustedHDmap.XLimits(2)+50]; 
% ylimits     = [adjustedHDmap.YLimits(1)-50 adjustedHDmap.YLimits(2)+50];
% zlimits     = [adjustedHDmap.ZLimits(1)-30 adjustedHDmap.ZLimits(2)+30];
ptCloudHD   = adjustedHDmap;
clearvars adjustedHDmap 
% player      = pcplayer(xlimits,ylimits,zlimits);
% title(player.Axes,'HD Point Cloud')
% xlabel(player.Axes,'X (m)');
% ylabel(player.Axes,'Y (m)');
% zlabel(player.Axes,'Z (m)');
% player.Axes.View = [0, 90];
% view(player, pcdownsample(ptCloudHD, 'gridAverage', 0.25)); 

%% Load Integrated Navigation Data: TC-INS/GNSS (IE)
initPOS                         = [23.0028474266, 120.214659453];
initPOS4WGS84                   = [initPOS(1), initPOS(2), 38.070];
[FileNameGT, PathNameGT, ~]     = uigetfile('*.txt', 'Please select your reference file (IE)');
f_ref                           = load([PathNameGT FileNameGT]);
t_gt                            = f_ref(1, 1);
LLF_NED_GT0                     = wgs2llf(f_ref, initPOS);
LLF_NED_GT                      = LLF_NED_GT0;
for i = 1:size(LLF_NED_GT,1)
    if (LLF_NED_GT(i,10)) < 0
        LLF_NED_GT(i,10) = LLF_NED_GT(i,10) + 360;
    else
        continue;
    end
end
clearvars LLF_NED_GT0

%% Output Configuration:
formatOut       = 'yymmdd';
headerFilename  = datestr(now,formatOut);
% insDir          = strcat(PathNameGT, headerFilename, '_INS_NAV', '.txt');
% ndtDir          = strcat(PathNameGT, headerFilename, '_NDT_NAV', '.txt');
% fid_ins         = fopen(insDir, 'w');
% fid_ndt         = fopen(ndtDir, 'w');
% outStr          = '%.8f\t%.20f\t%.20f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\t%.8f\n';

%% Main Process: 
currScan1       = 0; 
currScan2       = 0; 
currScan3       = 0; 
cntMatch33      = 0;
cntMatch31      = 0;
cntMatch32      = 0;
cntMatch312     = 0;
t_gt1           = t_gt;
t_gt2           = t_gt;
t_gt3           = t_gt;
start_gt_row1   = 1; 
start_gt_row2   = 1;
start_gt_row3   = 1;
cntTimeSlot     = 1;

while (hasFrame(veloReader3)) % HDL-64E (Middle or M-LiDAR)
    % Check an existing scan frame
    if ((currScan3+1) > veloReader3.NumberOfFrames)
        fprintf('End of Frame!');
        break;
    end
    currScan3   = currScan3 + 1;
    t_lidar3    = PCAP_UTC_TO_GPS_TIME(seconds(veloReader3.Timestamps(currScan3))*1000000, LiDAR_CFG3);
    fprintf('Primary LiDAR Frame No.%0.0f [%0.6f]\n', currScan3, t_lidar3);
    % Check an existing time of current LiDAR
    if (t_lidar3 < f_ref(1,1))
        continue;
    end	
	% Time Synchronization:
    while (t_gt3 <= t_lidar3)
        start_gt_row3   = start_gt_row3 + 1;
        t_gt3           = LLF_NED_GT(start_gt_row3, 1);        
    end
    % Interpolation:
    interWGSPose3                           = interpolatePOS(t_lidar3, f_ref, start_gt_row3);
    interLocalPose3                         = interpolatePOS(t_lidar3, LLF_NED_GT, start_gt_row3);
    tmpInterLocalpose3                      = TF_NED2ENU(interLocalPose3(:,1:10));
    LLF_ENU_GT.INT_SOL3(currScan3, 1:10)    = tmpInterLocalpose3;
    clearvars tmpInterLocalpose3 interLocalPose3
    % Get an initial navigation data:     
    nav_gt3      = getNavSolENU(LLF_ENU_GT.INT_SOL3, currScan3);
    r3_bn_gt     = nav_gt3.tran; 
    v3_bn_gt     = nav_gt3.vel; 
    a3_bn_gt     = nav_gt3.att';
    R3_n2b_gt    = nav_gt3.rotm;        
    R3_n2l_gt    = R3_n2b_gt * R3_b2l; 
    a3_ln_gt     = getEulerFromRotm(R3_n2l_gt, a3_bn_gt');
    R3_l2n_gt    = R3_n2l_gt';                             
    r3_ln_gt     = (R3_n2b_gt')*r3_lb + r3_bn_gt';          	

    if any(((t_lidar3+0.025) >= startTime) & (t_lidar3 <= endTime))
        cntMatch33 = cntMatch33 + 1;
        while (hasFrame(veloReader1)) % VLP-16 (Front or F-LiDAR) 
            % Check an existing scan frame
            if ((currScan1+1) > veloReader1.NumberOfFrames)
                fprintf('End of Frame!');
                return;
            end
            currScan1   = currScan1 + 1;
            t_lidar1    = PCAP_UTC_TO_GPS_TIME(seconds(veloReader1.Timestamps(currScan1))*1000000, LiDAR_CFG1);
            fprintf('Front LiDAR Frame No.%0.0f [%0.6f]\n', currScan1, t_lidar1);
            % Check an existing time of current LiDAR
            if (t_lidar1 < f_ref(1,1))
                continue;
            end	
            % Time Synchronization:
            while (t_gt1 <= t_lidar1)
                start_gt_row1   = start_gt_row1 + 1;
                t_gt1           = LLF_NED_GT(start_gt_row1, 1);        
            end
            % Interpolation:
            interLocalPose1                         = interpolatePOS(t_lidar1, LLF_NED_GT, start_gt_row1);
            tmpInterLocalpose1                      = TF_NED2ENU(interLocalPose1(:,1:10)); 
            LLF_ENU_GT.INT_SOL1(currScan1, 1:10)    = tmpInterLocalpose1;
            clearvars tmpInterLocalpose1 interLocalPose1
            % Get an initial navigation data:     
            nav_gt1      = getNavSolENU(LLF_ENU_GT.INT_SOL1, currScan1);
            r1_bn_gt     = nav_gt1.tran; 
            v1_bn_gt     = nav_gt1.vel; 
            a1_bn_gt     = nav_gt1.att';
            R1_n2b_gt    = nav_gt1.rotm;        
            R1_n2l_gt    = R1_n2b_gt * R1_b2l; 
            a1_ln_gt     = getEulerFromRotm(R1_n2l_gt, a1_bn_gt');
            R1_l2n_gt    = R1_n2l_gt';                             
            r1_ln_gt     = (R1_n2b_gt')*r1_lb + r1_bn_gt'; 

            % Synchronize current LiDAR-t1 w.r.t t3
            if (t_lidar1 > t_lidar3-0.03/lidarRate1)
                dt13 = t_lidar1 - t_lidar3;
                while (hasFrame(veloReader2)) % VLP-16 (Back or B-LiDAR) 
                    % Check an existing scan frame
                    if ((currScan2+1) > veloReader2.NumberOfFrames)
                        fprintf('End of Frame!');
                        return;
                    end
                    currScan2   = currScan2 + 1;
                    t_lidar2    = PCAP_UTC_TO_GPS_TIME(seconds(veloReader2.Timestamps(currScan2))*1000000, LiDAR_CFG2);
                    fprintf('Back LiDAR Frame No.%0.0f [%0.6f]\n', currScan2, t_lidar2);
                    % Check an existing time of current LiDAR
                    if (t_lidar2 < f_ref(1,1))
                        continue;
                    end	
                    % Time Synchronization:
                    while (t_gt2 <= t_lidar2)
                        start_gt_row2   = start_gt_row2 + 1;
                        t_gt2           = LLF_NED_GT(start_gt_row2, 1);        
                    end
                    % Interpolation:
                    interLocalPose2                         = interpolatePOS(t_lidar2, LLF_NED_GT, start_gt_row2);
                    tmpInterLocalpose2                      = TF_NED2ENU(interLocalPose2(:,1:10));
                    LLF_ENU_GT.INT_SOL2(currScan2, 1:10)    = tmpInterLocalpose2;
                    clearvars tmpInterLocalpose2 interLocalPose2
                    % Get an initial navigation data:     
                    nav_gt2      = getNavSolENU(LLF_ENU_GT.INT_SOL2, currScan2);
                    r2_bn_gt     = nav_gt2.tran; 
                    v2_bn_gt     = nav_gt2.vel; 
                    a2_bn_gt     = nav_gt2.att';
                    R2_n2b_gt    = nav_gt2.rotm;        
                    R2_n2l_gt    = R2_n2b_gt * R2_b2l; 
                    a2_ln_gt     = getEulerFromRotm(R2_n2l_gt, a2_bn_gt');
                    R2_l2n_gt    = R2_n2l_gt';                             
                    r2_ln_gt     = (R2_n2b_gt')*r2_lb + r2_bn_gt'; 
                    % Synchronize current LiDAR-t2 w.r.t t3
                    if (t_lidar2 > t_lidar3-0.03/lidarRate2)
                        dt23        = t_lidar2 - t_lidar3;
                        matchTime   = [currScan3, currScan1, currScan2; ...
                                       t_lidar3, t_lidar1, t_lidar2; ...
                                       (t_lidar3-t_lidar3), (t_lidar1-t_lidar3), (t_lidar2-t_lidar3)];
                        fprintf('\n -----------------LiDARs-Time Synchronization-----------------');
                        fprintf('\n                  Prime-LiDAR       1st-LiDAR        2nd-LiDAR');    
                        fprintf('\n    Scan No.       %0.3f        %0.3f        %0.3f', matchTime(1,1), matchTime(1,2), matchTime(1,3));
                        fprintf('\n    Timestamp     %0.3f       %0.3f       %0.3f', matchTime(2,1), matchTime(2,2), matchTime(2,3));
                        fprintf('\n    Timelapse        %0.3f            %0.3f            %0.3f', matchTime(3,1), matchTime(3,2), matchTime(3,3));
                        fprintf('\n ------------------------------------------------------------\n');

                        if (norm(matchTime(3,:)) >= 0.01)
                            break;
                        end
                        cntMatch312  = cntMatch312 + 1;

                        % Load Raw Point Cloud Object:
                        clearvars ptCloudRaw1 ptCloudRaw2 ptCloudRaw3 ptCloudOrg1 ptCloudOrg2 ptCloudOrg3
                        ptCloudRaw1  = readFrame(veloReader1,currScan1);
                        ptCloudRaw2  = readFrame(veloReader2,currScan2);
                        ptCloudRaw3  = readFrame(veloReader3,currScan3);

                        % Organize Raw Point Clouds:
                        ptCloudOrg1  = ptCloudOrgConversion(ptCloudRaw1, LiDAR_CFG1);
                        ptCloudOrg2  = ptCloudOrgConversion(ptCloudRaw2, LiDAR_CFG2);
                        ptCloudOrg3  = ptCloudOrgConversion(ptCloudRaw3, LiDAR_CFG3);

                        % Direct Geo-reference Point Cloud:
                        ptCloudDG3          = pcalign(ptCloudOrg3, rigid3d(R3_n2l_gt, r3_ln_gt'));
                        ROI3                = [ptCloudDG3.XLimits(1), ptCloudDG3.XLimits(2), ...
                                               ptCloudDG3.YLimits(1), ptCloudDG3.YLimits(2), ...
                                               ptCloudHD.ZLimits(1), ptCloudHD.ZLimits(2)];
                        ind3                = findPointsInROI(ptCloudHD, ROI3);
                        extractedHDmap      = select(ptCloudHD, ind3);
                        clearvars  ind3 ROI3 ptCloudDG3

                        % LiDAR Matching: M-LiDAR
                        movingScan3 = pcalign(scopePointCloud(ptCloudOrg3, minRange, maxRange3), rigid3d(R3_n2l_gt, [0,0,0]));
                        initTform3  = rigid3d(eye(3), r3_ln_gt');
                        [tform3, ~, ndtStatInfo3] = ndt_scan_matcher(movingScan3, extractedHDmap, LiDAR_CFG3.VoxelSizeNDT, ...
                                                    'InitialTransform', initTform3, 'MaxIterations', 50, ...
                                                    'Tolerance',  [0.01, 0.1]);
                        R3_n2l_ndt       = double((R3_n2l_gt)*(double(tform3.Rotation)));
                        a3_ln_ndt        = getEulerFromRotm(R3_n2l_ndt, a3_ln_gt);
                        R3_n2b_ndt       = double(R3_n2l_ndt)*double((R3_b2l'));
                        a3_bn_ndt        = getEulerFromRotm(R3_n2b_ndt, a3_bn_gt');
                        R3_b2l_ndt       = double(R3_n2l_ndt)/double(R3_n2b_gt);
                        a3_lb_ndt        = getEulerFromRotm(R3_b2l_ndt, [0;0;0]);
                        d_trans3         = double(tform3.Translation) - double(r3_ln_gt');
                        d_angle3         = (double(a3_ln_ndt')) - (double(a3_ln_gt'));
                        fprintf('\n-------------------- NDT Matching Performance --------------------\n');
                        fprintf('Diff. Translation (m.): %0.5f (E), %0.5f (N), %0.5f (U)\n', d_trans3);
                        fprintf('Diff. Rotation  (deg.): %0.5f (P), %0.5f (R), %0.5f (H)\n', d_angle3);
                        fprintf('------------------------------------------------------------------\n');
                        % Visualize the matched point clouds: Thresholding by positioning differences
                        if (cntMatch312 == 1) || (norm(d_trans3) >= 0.20) && (mod(cntMatch312,lidarRate3) == 0)
                            myVisualizePtCloud(currScan3, pcdownsample(extractedHDmap, 'gridAverage', 0.10), ...
                                pcalign(movingScan3, rigid3d(eye(3), r3_ln_gt')), pcalign(movingScan3, tform3), 'r', 'g');                 
                        end

                        % LiDAR Matching: F-LiDAR
                        movingScan1 = pcalign(scopePointCloud(ptCloudOrg1, minRange, maxRange1), rigid3d(R1_n2l_gt, [0,0,0]));
                        initTform1  = rigid3d(eye(3), r1_ln_gt');
                        [tform1, ~, ndtStatInfo1] = ndt_scan_matcher(movingScan1, extractedHDmap, LiDAR_CFG1.VoxelSizeNDT, ...
                                                    'InitialTransform', initTform1, 'MaxIterations', 50, ...
                                                    'Tolerance',  [0.01, 0.1]);
                        R1_n2l_ndt       = double((R1_n2l_gt)*(double(tform1.Rotation)));
                        a1_ln_ndt        = getEulerFromRotm(R1_n2l_ndt, a1_ln_gt);
                        R1_n2b_ndt       = double(R1_n2l_ndt)*double((R1_b2l'));
                        a1_bn_ndt        = getEulerFromRotm(R1_n2b_ndt, a1_bn_gt');
                        R1_b2l_ndt       = double(R1_n2l_ndt)/double(R1_n2b_gt);
                        a1_lb_ndt        = getEulerFromRotm(R1_b2l_ndt, [0;0;0]);
                        d_trans1         = double(tform1.Translation) - double(r1_ln_gt');
                        d_angle1         = (double(a1_ln_ndt')) - (double(a1_ln_gt'));
                        fprintf('\n-------------------- NDT Matching Performance --------------------\n');
                        fprintf('Diff. Translation (m.): %0.5f (E), %0.5f (N), %0.5f (U)\n', d_trans1);
                        fprintf('Diff. Rotation  (deg.): %0.5f (P), %0.5f (R), %0.5f (H)\n', d_angle1);
                        fprintf('------------------------------------------------------------------\n');
                        % Visualize the matched point clouds: Thresholding by positioning differences
                        if (cntMatch312 == 1) || (norm(d_trans1) >= 0.20) && (mod(cntMatch312,lidarRate1) == 0)
                            myVisualizePtCloud(currScan1, pcdownsample(extractedHDmap, 'gridAverage', 0.10), ...
                                pcalign(movingScan1, rigid3d(eye(3), r1_ln_gt')), pcalign(movingScan1, tform1), 'r', 'g');                 
                        end

                        % LiDAR Matching: B-LiDAR
                        movingScan2 = pcalign(scopePointCloud(ptCloudOrg2, minRange, maxRange2), rigid3d(R2_n2l_gt, [0,0,0]));
                        initTform2  = rigid3d(eye(3), r2_ln_gt');
                        [tform2, ~, ndtStatInfo2] = ndt_scan_matcher(movingScan2, extractedHDmap, LiDAR_CFG2.VoxelSizeNDT, ...
                                                    'InitialTransform', initTform2, 'MaxIterations', 50, ...
                                                    'Tolerance',  [0.01, 0.1]);
                        R2_n2l_ndt       = double((R2_n2l_gt)*(double(tform2.Rotation)));
                        a2_ln_ndt        = getEulerFromRotm(R2_n2l_ndt, a2_ln_gt);
                        R2_n2b_ndt       = double(R2_n2l_ndt)*double((R2_b2l'));
                        a2_bn_ndt        = getEulerFromRotm(R2_n2b_ndt, a2_bn_gt');
                        R2_b2l_ndt       = double(R2_n2l_ndt)/double(R2_n2b_gt);
                        a2_lb_ndt        = getEulerFromRotm(R2_b2l_ndt, [0;0;0]);
                        d_trans2         = double(tform2.Translation) - double(r2_ln_gt');
                        d_angle2         = (double(a2_ln_ndt')) - (double(a2_ln_gt'));
                        fprintf('\n-------------------- NDT Matching Performance --------------------\n');
                        fprintf('Diff. Translation (m.): %0.5f (E), %0.5f (N), %0.5f (U)\n', d_trans2);
                        fprintf('Diff. Rotation  (deg.): %0.5f (P), %0.5f (R), %0.5f (H)\n', d_angle2);
                        fprintf('------------------------------------------------------------------\n');
                        % Visualize the matched point clouds: Thresholding by positioning differences
                        if (cntMatch312 == 1) || (norm(d_trans2) >= 0.20) && (mod(cntMatch312,lidarRate2) == 0)
                            myVisualizePtCloud(currScan2, pcdownsample(extractedHDmap, 'gridAverage', 0.10), ...
                                pcalign(movingScan2, rigid3d(eye(3), r2_ln_gt')), pcalign(movingScan2, tform2), 'r', 'g');                 
                        end

                        % LiDAR Matching: F-to-M-LiDAR
                        movingScan13 = pcalign(movingScan1, rigid3d(tform1.Rotation, [0,0,0]));
                        initTform13  = rigid3d(eye(3), tform1.Translation);
                        [tform13, ~, ndtStatInfo13] = ndt_scan_matcher(movingScan13, pcalign(movingScan3, tform3), ...
                                            LiDAR_CFG1.VoxelSizeNDT,'InitialTransform', initTform13, 'MaxIterations', 50, ...
                                            'Tolerance',  [0.01, 0.1]);
                        R_l3_l1_ndt       = double((initTform13.Rotation)*(double(tform13.Rotation)));
                        a_l1_l3_ndt       = getEulerFromRotm(R_l3_l1_ndt, [0;0;0]);
                        R1_n2l_ndt_upt    = double((R1_n2l_ndt)*(double(R_l3_l1_ndt)));
                        a1_ln_ndt_upt     = getEulerFromRotm(R1_n2l_ndt_upt, a1_ln_ndt);
                        R1_b2l1_ndt       = double((R1_b2l_ndt)*(double(R_l3_l1_ndt)));
                        a_l1_b_ndt        = getEulerFromRotm(R1_b2l1_ndt, [0;0;0]);
                        d_trans13         = double(tform13.Translation) - double(initTform13.Translation);
                        d_angle13         = (double(a_l1_l3_ndt')) - (double([0;0;0]'));
                        fprintf('\n-------------------- NDT Matching Performance --------------------\n');
                        fprintf('Diff. Translation (m.): %0.5f (E), %0.5f (N), %0.5f (U)\n', d_trans13);
                        fprintf('Diff. Rotation  (deg.): %0.5f (P), %0.5f (R), %0.5f (H)\n', d_angle13);
                        fprintf('------------------------------------------------------------------\n');
                        % Visualize the matched point clouds: Thresholding by positioning differences
                        if (cntMatch312 == 1) || (norm(d_trans13) >= 0.15) && (mod(cntMatch312,lidarRate1) == 0)
                            myVisualizePtCloud(currScan1, pcalign(movingScan3, tform3), ...
                                pcalign(movingScan13, initTform13), pcalign(movingScan13, tform13), 'r', 'g');                 
                        end

                        % LiDAR Matching: B-to-M-LiDAR
                        movingScan23 = pcalign(movingScan2, rigid3d(tform2.Rotation, [0,0,0]));
                        initTform23  = rigid3d(eye(3), tform2.Translation);
                        [tform23, ~, ndtStatInfo23] = ndt_scan_matcher(movingScan23, pcalign(movingScan3, tform3), ...
                                            LiDAR_CFG2.VoxelSizeNDT, 'InitialTransform', initTform23, 'MaxIterations', 50, ...
                                            'Tolerance',  [0.01, 0.1]);
                        R_l3_l2_ndt       = double((initTform23.Rotation)*(double(tform23.Rotation)));
                        a_l2_l3_ndt       = getEulerFromRotm(R_l3_l2_ndt, [0;0;0]);
                        R2_n2l_ndt_upt    = double((R2_n2l_ndt)*(double(R_l3_l2_ndt)));
                        a2_ln_ndt_upt     = getEulerFromRotm(R2_n2l_ndt_upt, a2_ln_ndt);
                        R2_b2l2_ndt       = double((R2_b2l_ndt)*(double(R_l3_l2_ndt)));
                        a_l2_b_ndt        = getEulerFromRotm(R2_b2l2_ndt, [0;0;0]);
                        d_trans23         = double(tform23.Translation) - double(initTform23.Translation);
                        d_angle23         = (double(a_l2_l3_ndt')) - (double([0;0;0]'));
                        fprintf('\n-------------------- NDT Matching Performance --------------------\n');
                        fprintf('Diff. Translation (m.): %0.5f (E), %0.5f (N), %0.5f (U)\n', d_trans23);
                        fprintf('Diff. Rotation  (deg.): %0.5f (P), %0.5f (R), %0.5f (H)\n', d_angle23);
                        fprintf('------------------------------------------------------------------\n');
                        % Visualize the matched point clouds: Thresholding by positioning differences
                        if (cntMatch312 == 1) || (norm(d_trans23) >= 0.15) && (mod(cntMatch312,lidarRate2) == 0)
                            myVisualizePtCloud(currScan2, pcalign(movingScan3, tform3), ...
                                pcalign(movingScan23, initTform23), pcalign(movingScan23, tform23), 'r', 'g');                 
                        end
    
                        % Fault detection & isolation:
                        if (cntMatch312 > 1)
                            % Translation offset:
                            del_r1_ln = double(tform13.Translation) - prev_r1_ln;
                            del_r2_ln = double(tform23.Translation) - prev_r2_ln;
                            del_r3_ln = double(tform3.Translation) - prev_r3_ln;
                            del_r1_bn = double(r1_bn_gt) - prev_r1_bn;
                            del_r2_bn = double(r2_bn_gt) - prev_r2_bn;
                            del_r3_bn = double(r3_bn_gt) - prev_r3_bn;

                            % Rotation offset:
                            del_a1_ln = a1_ln_ndt_upt' - prev_a1_ln;
                            del_a2_ln = a2_ln_ndt_upt' - prev_a2_ln;
                            del_a3_ln = a3_ln_ndt' - prev_a3_ln;
                            del_a1_bn = double(a1_bn_gt) - prev_a1_bn;
                            del_a2_bn = double(a2_bn_gt) - prev_a2_bn;
                            del_a3_bn = double(a3_bn_gt) - prev_a3_bn;

                            % Lever-arms & bore-sight angles:
                            del_r1_lb = (double(R3_n2b_gt) * (double(tform13.Translation-r1_bn_gt)'))' - prev_r1_lb;
                            del_r2_lb = (double(R3_n2b_gt) * (double(tform23.Translation-r2_bn_gt)'))' - prev_r2_lb;
                            del_r3_lb = (double(R3_n2b_gt) * (double(tform3.Translation-r3_bn_gt)'))' - prev_r3_lb;
                            del_a1_lb = double(a_l1_b_ndt') - double(prev_a1_lb');
                            del_a2_lb = double(a_l2_b_ndt') - double(prev_a2_lb');
                            del_a3_lb = double(a3_lb_ndt') - double(prev_a3_lb');

                            if (abs(norm(del_r3_ln)-norm(del_r3_bn)) >= 0.05) || (abs(norm(del_a3_ln)-norm(del_a3_bn)) >= 0.5) ...
                                    || (abs(norm(del_r3_lb)) >= 0.05) || (abs(norm(del_a3_lb)) >= 0.5)
                                cntMatch312 = cntMatch312 - 1;
                                fprintf('\nOutlier FOUND!\n');
                                break;
                            end
                            if (abs(norm(del_r1_ln)-norm(del_r1_bn)) >= 0.05) || (abs(norm(del_a1_ln)-norm(del_a1_bn)) >= 0.5) ...
                                    || (abs(norm(del_r1_lb)) >= 0.05) || (abs(norm(del_a1_lb)) >= 0.5)
                                cntMatch312 = cntMatch312 - 1;
                                fprintf('\nOutlier FOUND!\n');
                                break;
                            end
                            if (abs(norm(del_r2_ln)-norm(del_r2_bn)) >= 0.05) || (abs(norm(del_a2_ln)-norm(del_a2_bn)) >= 0.5) ...
                                    || (abs(norm(del_r2_lb)) >= 0.05) || (abs(norm(del_a2_lb)) >= 0.5)
                                cntMatch312 = cntMatch312 - 1;
                                fprintf('\nOutlier FOUND!\n');
                                break;
                            end
                        end

                        if (norm(d_trans3) <= 0.15) && (abs(d_angle3(3)) <= 2.5)
                            % Store IMU's translation vector & rotation matrix: 
                            imuInfo.r1_bn(cntMatch312,1:3)      = double(r1_bn_gt);
                            imuInfo.r2_bn(cntMatch312,1:3)      = double(r2_bn_gt);
                            imuInfo.r3_bn(cntMatch312,1:3)      = double(r3_bn_gt);
                            imuInfo.a1_bn(cntMatch312,1:3)      = double(a1_bn_gt);
                            imuInfo.a2_bn(cntMatch312,1:3)      = double(a2_bn_gt);
                            imuInfo.a3_bn(cntMatch312,1:3)      = double(a3_bn_gt);
                            imuInfo.R1_n2b{cntMatch312,1}       = double(R1_n2b_gt);
                            imuInfo.R2_n2b{cntMatch312,1}       = double(R2_n2b_gt);
                            imuInfo.R3_n2b{cntMatch312,1}       = double(R3_n2b_gt);
                            % Store LiDAR's translation vector & rotation matrix:
                            lidarInfo.r1_ln(cntMatch312,1:3)    = double(tform13.Translation);
                            lidarInfo.r2_ln(cntMatch312,1:3)    = double(tform23.Translation);
                            lidarInfo.r3_ln(cntMatch312,1:3)    = double(tform3.Translation);
                            lidarInfo.a1_ln(cntMatch312,1:3)    = a1_ln_ndt_upt';
                            lidarInfo.a2_ln(cntMatch312,1:3)    = a2_ln_ndt_upt';
                            lidarInfo.a3_ln(cntMatch312,1:3)    = a3_ln_ndt';
                            lidarInfo.R1_n2l{cntMatch312,1}     = R1_n2l_ndt_upt;
                            lidarInfo.R2_n2l{cntMatch312,1}     = R2_n2l_ndt_upt;
                            lidarInfo.R3_n2l{cntMatch312,1}     = double(R3_n2l_ndt);
                            % Store LiDAR mounting parameters:
                            lidarInfo.r1_lb(cntMatch312,1:3)    = (double(R3_n2b_gt) * (double(tform13.Translation-r1_bn_gt)'))';
                            lidarInfo.r2_lb(cntMatch312,1:3)    = (double(R3_n2b_gt) * (double(tform23.Translation-r2_bn_gt)'))';
                            lidarInfo.r3_lb(cntMatch312,1:3)    = (double(R3_n2b_gt) * (double(tform3.Translation-r3_bn_gt)'))';
                            lidarInfo.a1_lb(cntMatch312,1:3)    = a_l1_b_ndt';
                            lidarInfo.a2_lb(cntMatch312,1:3)    = a_l2_b_ndt';
                            lidarInfo.a3_lb(cntMatch312,1:3)    = a3_lb_ndt';
                            lidarInfo.R1_b2l{cntMatch312,1}     = R1_b2l1_ndt;
                            lidarInfo.R2_b2l{cntMatch312,1}     = R2_b2l2_ndt;
                            lidarInfo.R3_b2l{cntMatch312,1}     = R3_b2l_ndt;

                            la_before   = [r1_lb', r2_lb', r3_lb'];
                            bs_before   = [bs10', bs20', bs30'];
                            la_after    = [(double(R3_n2b_gt) * (double(tform13.Translation-r1_bn_gt)'))', ...
                                            (double(R3_n2b_gt) * (double(tform23.Translation-r2_bn_gt)'))', ...
                                            (double(R3_n2b_gt) * (double(tform3.Translation-r3_bn_gt)'))'];
                            bs_after    = [a_l1_b_ndt', a_l2_b_ndt', a3_lb_ndt'];
                            fprintf('\n-------------------------------------- LiDAR-Mounting Parameters  --------------------------------------');
                            fprintf('\n----------------------------------------- (Before Calibration) -----------------------------------------');
                            fprintf('\n    [X_RIGHT; Y_FWD; Z_UP]        1st LiDAR                 2nd LiDAR             3rd LiDAR');    
                            fprintf('\n          Lever-arms (m.):  [%0.3f; %0.3f; %0.3f]   [%0.3f; %0.3f; %0.3f]   [%0.3f; %0.3f; %0.3f]', la_before);
                            fprintf('\n Bore-sight angles (deg.):  [%0.3f; %0.3f; %0.3f]   [%0.3f; %0.3f; %0.3f]   [%0.3f; %0.3f; %0.3f]', bs_before);
                            fprintf('\n------------------------------------------ (After Calibration) -----------------------------------------');
                            fprintf('\n    [X_RIGHT; Y_FWD; Z_UP]        1st LiDAR               2nd LiDAR               3rd LiDAR');    
                            fprintf('\n          Lever-arms (m.):  [%0.3f; %0.3f; %0.3f]   [%0.3f; %0.3f; %0.3f]   [%0.3f; %0.3f; %0.3f]', la_after);
                            fprintf('\n Bore-sight angles (deg.):  [%0.3f; %0.3f; %0.3f]   [%0.3f; %0.3f; %0.3f]   [%0.3f; %0.3f; %0.3f]', bs_after);
                            fprintf('\n--------------------------------------------------------------------------------------------------------\n');

                            % For next epoch time:
                            % == Translation:
                            prev_r1_ln = double(tform13.Translation);
                            prev_r2_ln = double(tform23.Translation);
                            prev_r3_ln = double(tform3.Translation);
                            prev_r1_bn = double(r1_bn_gt);
                            prev_r2_bn = double(r2_bn_gt);
                            prev_r3_bn = double(r3_bn_gt);
                            % == Rotation:
                            prev_a1_ln = double(a1_ln_ndt_upt');
                            prev_a2_ln = double(a2_ln_ndt_upt');
                            prev_a3_ln = double(a3_ln_ndt');
                            prev_a1_bn = double(a1_bn_gt);
                            prev_a2_bn = double(a2_bn_gt);
                            prev_a3_bn = double(a3_bn_gt);
                            % == Lever-arms & bore-sight angles:
                            prev_r1_lb = lidarInfo.r1_lb(cntMatch312,1:3);
                            prev_r2_lb = lidarInfo.r2_lb(cntMatch312,1:3);
                            prev_r3_lb = lidarInfo.r3_lb(cntMatch312,1:3);
                            prev_a1_lb = a_l1_b_ndt;
                            prev_a2_lb = a_l2_b_ndt;
                            prev_a3_lb = a3_lb_ndt;
                        else
                            cntMatch312 = cntMatch312 - 1;
                            fprintf('\nOutlier FOUND!\n');
                            break;
                        end
                        break;
                    else
                        continue;
                    end
                end
                break;
            else
                continue;
            end
        end
    else
       continue;
    end
end

% Save navigation data for lidar calibration:
navInfo.imu     = imuInfo;
navInfo.lidar   = lidarInfo;
[fname, pname]  = uiputfile('*.mat', 'Save Nav. Info.');
[~,name,ext]    = fileparts(fname);
save(strcat(pname, name, ext),'-v7.3', 'navInfo');




















