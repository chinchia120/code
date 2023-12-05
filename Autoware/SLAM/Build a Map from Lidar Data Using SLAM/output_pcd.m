%% =============== Function of output_pcd =============== %%
function output_pcd(ptCloudMap, FileName_output)
    % Assuming ptCloud is your point cloud object and 'intensity' is a vector of intensity values
    minI = min(double(ptCloudMap.Intensity));
    maxI = max(double(ptCloudMap.Intensity));
    
    % Normalize intensity to 0-1
    normI = (double(ptCloudMap.Intensity)-minI) / (maxI-minI);
    
    % Map intensity to an RGB color (here we use a grayscale colormap, but others can be used)
    colorMap = uint8(255*repmat(normI, 1, 3));
    
    % Assign color to point cloud
    ptCloudMap.Color = colorMap;
    
    % Write to PLY file
    pcwrite(ptCloudMap, FileName_output);
end

