% Coordinate Transformation from WGS84 to TWD97
function out_file = wgs84Totwd97(latlon_data)
out_file = latlon_data;

% Get parameters of WGS84-to-TWD97:
params = get_wgs84_to_twd97_params();

% Transformation: WGS84 ---> TWD97
[row, col] = size(latlon_data);
if ~isequal(col, 2)
    return;
end
for i = 1:row
    lat = deg2rad(latlon_data(i,1));
    lon = deg2rad(latlon_data(i,2));
    [x, y] = transform_wgs84_to_twd97(lat, lon, params);

    out_file(i, 1) = x;
    out_file(i, 2) = y;
end
end