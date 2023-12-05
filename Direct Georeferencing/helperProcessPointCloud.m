%% =============== Function of helperProcessPointCloud =============== %%
function ptCloud = helperProcessPointCloud(ptCloudIn, method)
    arguments
        ptCloudIn (1,1) pointCloud
        method string {mustBeMember(method, ["planefit","rangefloodfill"])} = "rangefloodfill"
    end
    
    isOrganized = ~ismatrix(ptCloudIn.Location);
    
    if (method=="rangefloodfill" && isOrganized) 
        % Segment ground using floodfill on range image
        groundFixedIdx = segmentGroundFromLidarData(ptCloudIn, "ElevationAngleDelta", 11);
    else
        % Segment ground as the dominant plane with reference normal
        % vector pointing in positive z-direction
        maxDistance = 0.4;
        maxAngularDistance = 5;
        referenceVector = [0, 0, 1];
    
        [~, groundFixedIdx] = pcfitplane(ptCloudIn, maxDistance, referenceVector, maxAngularDistance);
    end
    
    if isOrganized
        groundFixed = false(size(ptCloudIn.Location, 1), size(ptCloudIn.Location, 2));
    else
        groundFixed = false(ptCloudIn.Count, 1);
    end
    groundFixed(groundFixedIdx) = true;
    
    % Segment ego vehicle as points within a given radius of sensor
    sensorLocation = [0, 0, 0];
    radius = 4;
    egoFixedIdx = findNeighborsInRadius(ptCloudIn, sensorLocation, radius);
    
    if isOrganized
        egoFixed = false(size(ptCloudIn.Location, 1), size(ptCloudIn.Location, 2));
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