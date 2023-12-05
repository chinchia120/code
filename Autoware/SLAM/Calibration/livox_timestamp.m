function [t_week,t_lidar] = livox_timestamp(livox_pcds,idx)
    t_lidar = string(livox_pcds.Files{idx});
    t_lidar = split(t_lidar,'.');
    t_lidar = split(t_lidar(1),'\');
    t_lidar = str2double(t_lidar(end))*10^-9;
    utc = unix2time(t_lidar);
    gpst = utc2gpst(utc);
    [t_week, t_lidar] = time2gpst(gpst);
end

