%% ========== Setup ========== %%
clc;
clear all;
close all;

%% ========== Initial Value ========== %%
student_id = 'P66134111';
longitude = mod(str2num(extractAfter(student_id, 1)), 180);             % deg
latitude = mod(str2num(extractBefore(reverse(student_id), 9)), 90);     % deg
altitude = 0;                                                           % m
recPos = [latitude longitude altitude];                                 % [deg deg m]
maskAngle = 0;                                                          % deg
gnssFileType = "RINEX";                                                 % s

%% ========== Display the GNSS Receiver Position ========== %%
geoplot(latitude, longitude, "x");
title("Receiver Position");

% ========== Get Satellite Orbital Parameters ========== %%
satSys = "GPS";
[navmsg, satIDs] = exampleHelperParseGNSSFile(gnssFileType, satSys);

%% ========== Skyplot of GPS ========== %%
startTime = datetime(2021,06,24,04,00,00);
secondsPerHour = 3600;
numHours = 24;
dt = 60*60;                 
day = 365*2+281;
timeElapsed = 0: dt :(secondsPerHour*numHours)*day+(secondsPerHour*20);
t = startTime+seconds(timeElapsed);

sp = skyplot([], [], MaskElevation=maskAngle);

for ii = 1: numel(t)
    satPos = gnssconstellation(t(ii), navmsg, GNSSFileType=gnssFileType);
    [az(ii,:), el(ii,:), vis(ii,:)] = lookangles(recPos, satPos, maskAngle);
    set(sp, AzimuthData=az(ii, vis(ii,:)), ...
            ElevationData=el(ii, vis(ii,:)), ...
            LabelData=satIDs(vis(ii,:)));
    title(sprintf('Skyplot of GPS at %s', t(ii)));
    drawnow limitrate;
end

% Save Plot
saveas(gcf, sprintf('2024-04-01-000000-1.fig'));

%% ========== Skyplot of Satellite ========== %%
allSatSys = ["GPS","Galileo","GLONASS","BeiDou","NavIC","QZSS","SBAS"];
satLetter = ["G","E","R","C","I","J","S"];

sp = skyplot([],[],MaskElevation=maskAngle);

for ii = 1:numel(allSatSys)
    satSys = allSatSys(ii);

    % Get the satellite positions of the current satellite system
    [navmsg,satIDs] = exampleHelperParseGNSSFile("RINEX",satSys);
    % if(satSys == "GLONASS" || satSys == "BeiDou" || satSys == "NavIC" || satSys == "SBAS"); continue; end
    satPos = gnssconstellation(startTime,navmsg);
    [az,el,vis] = lookangles(recPos,satPos,maskAngle);
    
    % Combine the satellite system symbol letter with the satellite ID
    satIDLabel = arrayfun(@(x) sprintf("%c%02d",satLetter(ii),x),satIDs);
    
    % Create a categorical array to associate the current values with the current satellite system in the skyplot
    satGroup = categorical(repmat(ii,numel(satIDLabel),1),1:numel(allSatSys),allSatSys);
    
    % Update the skyplot
    set(sp, AzimuthData=[sp.AzimuthData(:); az(vis)], ...
            ElevationData=[sp.ElevationData(:); el(vis)], ...
            LabelData=[sp.LabelData(:); satIDLabel(vis)], ...
            GroupData=[sp.GroupData; satGroup(vis)]);
    title(sprintf('Skyplot of GPS at %s', t(ii)));
end

% Add a legend to the skyplot
legend;

% Save Plot
saveas(gcf, sprintf('2024-04-01-000000-2.fig'));

%% ========== Helper Function ========== %%
function [navmsg,satIDs] = exampleHelperParseGNSSFile(gnssFileType,satSys)
    switch gnssFileType
        case "RINEX"
            switch satSys
                case "GPS"
                    file = "GODS00USA_R_20211750000_01D_GN.rnx";
                    navmsg = rinexread(file);
                    % For RINEX files, extract GPS data and use only the first entry for each satellite.
                    gpsData = navmsg.GPS;
                    [~,idx] = unique(gpsData.SatelliteID);
                    navmsg = gpsData(idx,:);
                case "Galileo"
                    file = "GODS00USA_R_20211750000_01D_EN.rnx";
                    navmsg = rinexread(file);
                    % For RINEX files, extract Galileo data and use only the first entry for each satellite.
                    galData = navmsg.Galileo;
                    [~,idx] = unique(galData.SatelliteID);
                    navmsg = galData(idx,:);
                case "GLONASS"
                    file = "GODS00USA_R_20211750000_01D_RN.rnx";
                    navmsg = rinexread(file);
                    % For RINEX files, extract GLONASS data and use only the first entry for each satellite.
                    gloData = navmsg.GLONASS;
                    [~,idx] = unique(gloData.SatelliteID);
                    navmsg = gloData(idx,:);
                case "BeiDou"
                    file = "GODS00USA_R_20211750000_01D_CN.rnx";
                    navmsg = rinexread(file);
                    % For RINEX files, extract BeiDou data and use only the first entry for each satellite.
                    bdsData = navmsg.BeiDou;
                    [~,idx] = unique(bdsData.SatelliteID);
                    navmsg = bdsData(idx,:);
                case "NavIC"
                    file = "ARHT00ATA_R_20211750000_01D_IN.rnx";
                    navmsg = rinexread(file);
                    % For RINEX files, extract NavIC data and use only the first entry for each satellite.
                    irnData = navmsg.NavIC;
                    [~,idx] = unique(irnData.SatelliteID);
                    navmsg = irnData(idx,:);
                case "QZSS"
                    file = "ARHT00ATA_R_20211750000_01D_JN.rnx";
                    navmsg = rinexread(file);
                    % For RINEX files, extract QZSS data and use only the first entry for each satellite.
                    qzsData = navmsg.QZSS;
                    [~,idx] = unique(qzsData.SatelliteID);
                    navmsg = qzsData(idx,:);
                case "SBAS"
                    file = "GOP600CZE_R_20211750000_01D_SN.rnx";
                    navmsg = rinexread(file);
                    % For RINEX files, extract SBAS data and use only the first entry for each satellite.
                    sbaData = navmsg.SBAS;
                    [~,idx] = unique(sbaData.SatelliteID);
                    navmsg = sbaData(idx,:);
            end
            satIDs = navmsg.SatelliteID;
        case "SEM"
            file = "semalmanac_2021-6-22.al3";
            navmsg = semread(file);
            satIDs = navmsg.PRNNumber;
        case "YUMA"
            file = "yumaalmanac_2021-6-22.alm";
            navmsg = yumaread(file);
            satIDs = navmsg.PRN;
    end
end