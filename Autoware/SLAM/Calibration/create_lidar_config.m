%==========================================================================
% Purposes: Create LiDAR-configuration from an input LiDAR file
% Function: create_lidar_config
% Developed & Coded by: Surachet Srinara, PhD candidate
% Date: November 8, 2021 - 03:40PM
% Copyright@Chet2021
%========================================================================== 
function [LiDAR_CFG1, veloReader1] = create_lidar_config(pcPathName1, pcFileName1)
splitStr1 = strsplit(pcFileName1, {'-','_','.'});
laserModel1 = strcat(splitStr1{8}, splitStr1{9});
veloReader1 = velodyneFileReader([pcPathName1 pcFileName1],laserModel1);
LiDAR_CFG1.numFrames = veloReader1.NumberOfFrames;

dateTime1 = zeros(1,6);
for i = 1:6
    dateTime1(i) = str2double(splitStr1{i});
end
t1 = datetime(dateTime1(1),dateTime1(2),dateTime1(3));
d1 = day(t1,'dayofweek');
LiDAR_CFG1.DOW = d1-1;                                % Day of week (DOW)
LiDAR_CFG1.HOUR = dateTime1(4);                       % Hours
LiDAR_CFG1.MINS = dateTime1(5);                       % Minutes

TimeZone = 'Asia/Taipei';
T = timezones('Asia');
for k = 1:size(T,1)
    if isequal(T.Name{k}, TimeZone)
        your_current_timezone = T.UTCOffset(k);
        break;
    else
        continue;
    end
end

% Default setting:
LiDAR_CFG1.veloReader = veloReader1;
LiDAR_CFG1.TimeZone_NOW = your_current_timezone;
LiDAR_CFG1.GPS_LEAP_SEC = 18;                           % GPS Leap Seconds
LiDAR_CFG1.SEC_IN_ONE_DAY = 24*60*60;                   % Seconds in one day
LiDAR_CFG1.SEC_IN_ONE_HOUR = 60*60;                     % Seconds in one hour
LiDAR_CFG1.UTC_DELAY = 0.00;                            % UTC delay
LiDAR_CFG1.lidar_min_range = 1.5;                       % Minimum range
LiDAR_CFG1.lidar_max_range = 70;                        % Maximum range
LiDAR_CFG1.VoxelSizeNDT = 3;                            % NDT voxel size (m.)
LiDAR_CFG1.FilterSizeNDT = 0.10;                        % NDT filter size (m.)

end