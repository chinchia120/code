%% =============== Function of helperProcessPointCloud_SetRange =============== %%
function ptCloud = helperProcessPointCloud_SetRange(ptCloudIn, method)
    % Initial Setup
    arguments
        ptCloudIn (1,1) pointCloud
        method string {mustBeMember(method, ["planefit", "rangefloodfill"])} = "rangefloodfill"
    end
    
    isOrganized = ~ismatrix(ptCloudIn.Location);
    
    if (method == "rangefloodfill" && isOrganized) 
    % Segment Ground Using Floodfill on Range Image
        groundFixedIdx = segmentGroundFromLidarData(ptCloudIn, "ElevationAngleDelta", 11);
    else
    % Segment Ground as the Dominant Plane with Reference Normal Vector Pointing in Positive Z-Direction
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
    
    % Segment Ego Vehicle as Points within a Given Radius of Sensor
    sensorLocation = [0, 0, 0];
    radius = 4;
    egoFixedIdx = findNeighborsInRadius(ptCloudIn, sensorLocation, radius);
    
    if isOrganized
        egoFixed = false(size(ptCloudIn.Location, 1), size(ptCloudIn.Location, 2));
    else
        egoFixed = false(ptCloudIn.Count, 1);
    end
    egoFixed(egoFixedIdx) = true;
    
    % Retain Subset of Point Cloud without Ground and Ego Vehicle
    if isOrganized
        indices = ~groundFixed & ~egoFixed;
    else
        indices = find(~groundFixed & ~egoFixed);
    end
    
    ptCloud_rm = pcdenoise(select(ptCloudIn, indices)); 

    % Select Range of Point Cloud
    roi = [-15, 15, -15, 15, min(ptCloud_rm.Location(:, 3)), max(ptCloud_rm.Location(:, 3))];
    indices = findPointsInROI(ptCloud_rm, roi);
    
    ptCloud = select(ptCloud_rm, indices);
end