% Unorganized to origanized conversion of point clouds using spherical projection
% Surachet Srinara
% May 22, 2022

function ptCloudOrg = ptCloudOrgConversion(ptCloudUnOrg, LiDAR_CFG)
% Specify the sensor parameters using lidarParameters function.
params = lidarParameters(LiDAR_CFG.vResolution, LiDAR_CFG.vFoV, ...
                        LiDAR_CFG.hResolution,"HorizontalFoV",LiDAR_CFG.hFoV);

% Convert the unorganized point cloud to organized format using the pcorganize function.
ptCloudOrg = pcorganize(ptCloudUnOrg, params);

% Display the intensity channel of the reconstructed organized point cloud. 
% Resize the image for better visualization.
% intensityChannel = ptCloudOrg.Intensity;
% intensityChannel = imresize(intensityChannel,'Scale',[3 1]);
% figure
% imshow(intensityChannel);
% xChannel = ptCloudOrg.Location(:,:,1);
% xChannel = imresize(xChannel,'Scale',[3 1]);
% figure
% imshow(xChannel);
% yChannel = ptCloudOrg.Location(:,:,2);
% yChannel = imresize(yChannel,'Scale',[3 1]);
% figure
% imshow(yChannel);
% zChannel = ptCloudOrg.Location(:,:,3);
% zChannel = imresize(zChannel,'Scale',[3 1]);
% figure
% imshow(zChannel);

% Display both the original organized point cloud and the reconstructed organized point cloud using 
% the helperShowUnorgAndOrgPair helper function, attached to this example as a supporting file.
% zoomFactor  = 2.5;
% display3    = myhelperShowUnorgAndOrgPair();
% display3.plotLidarScan(ptCloudUnOrg,ptCloudOrg,zoomFactor);

end