%% ========== Setup ========== %%
clc;
clear all;
close all;

%% ========== Specify Simulation Parameters ========== %%
gnssFileType = "RINEX";
startTime = datetime(2021,06,24,04,00,00);
numHours = 24;
dt =60; % s
latitude =42.3013162; % deg
longitude = -71.3782972; % deg
altitude = 50; % m
recPos = [latitude longitude altitude]; % [deg deg m]
maskAngle = 5; % deg

geoplot(latitude,longitude,"x")
title("Receiver Position")

satSys = "GPS";
if gnssFileType == "RINEX"
    satSys = "GPS";
end

%% ========== Get Satellite Orbital Parameters ========== %%
[navmsg,satIDs] = exampleHelperParseGNSSFile(gnssFileType,satSys);

%% ========== Generate Satellite Visibilities ========== %%
satPos = gnssconstellation(startTime,navmsg,GNSSFileType=gnssFileType);
[az,el,vis] = lookangles(recPos,satPos,maskAngle);
figure
skyplot(az(vis),el(vis),satIDs(vis),MaskElevation=maskAngle);

secondsPerHour = 3600;
timeElapsed = 0:dt:(secondsPerHour*numHours);
t = startTime+seconds(timeElapsed);

numSats = numel(satIDs);
numSamples = numel(t);
az = zeros(numSamples,numSats);
el = zeros(numSamples,numSats);
vis = false(numSamples,numSats);

sp = skyplot([],[],MaskElevation=maskAngle);

for ii = 1:numel(t)
    satPos = gnssconstellation(t(ii),navmsg,GNSSFileType=gnssFileType);
    [az(ii,:),el(ii,:),vis(ii,:)] = lookangles(recPos,satPos,maskAngle);
    set(sp,AzimuthData=az(ii,vis(ii,:)), ...
        ElevationData=el(ii,vis(ii,:)), ...
        LabelData=satIDs(vis(ii,:)))
    drawnow limitrate
end

%% ========== Plot Results ========== %%
visPlotData = double(vis);
visPlotData(visPlotData == false) = NaN; % Hide invisible satellites.
visPlotData = visPlotData + (0:numSats-1); % Add space to satellites to be stacked.
colors = colororder;

figure
plot(t,visPlotData,".",Color=colors(1,:))
yticks(1:numSats)
yticklabels(string(satIDs))
grid on
ylabel("Satellite ID")
xlabel("Time")
title("Satellite Visibility Chart")
axis tight

numVis = sum(vis,2);
figure
area(t,numVis)
grid on
xlabel("Time")
ylabel("Number of satellites visible")
title("Number of GNSS satellites visible")
axis tight

%% ========== Multi-Constellation Visibilities ========== %%
allSatSys = ["GPS","Galileo","GLONASS","BeiDou","NavIC","QZSS","SBAS"];
satLetter = ["G","E","R","C","I","J","S"];

sp = skyplot([],[],MaskElevation=maskAngle);

for ii = 1:numel(allSatSys)
    satSys = allSatSys(ii);

    % Get the satellite positions of the current satellite system.
    [navmsg,satIDs] = exampleHelperParseGNSSFile("RINEX",satSys);
    if(satSys == "GLONASS" || satSys == "BeiDou" || satSys == "NavIC" || satSys == "SBAS"); continue; end
    satPos = gnssconstellation(startTime,navmsg);
    [az,el,vis] = lookangles(recPos,satPos,maskAngle);
    
    % Combine the satellite system symbol letter with the satellite ID.
    satIDLabel = arrayfun(@(x) sprintf("%c%02d",satLetter(ii),x),satIDs);
    
    % Create a categorical array to associate the current values with the current satellite system in the skyplot.
    satGroup = categorical(repmat(ii,numel(satIDLabel),1),1:numel(allSatSys),allSatSys);
    
    % Update the skyplot.
    set(sp, ...
        AzimuthData=[sp.AzimuthData(:); az(vis)], ...
        ElevationData=[sp.ElevationData(:); el(vis)], ...
        LabelData=[sp.LabelData(:); satIDLabel(vis)], ...
        GroupData=[sp.GroupData; satGroup(vis)]);
    %break;
end

% Add a legend to the skyplot.
legend

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