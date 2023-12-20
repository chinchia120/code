% Get WGS84-to-TWD97 parameters
function cfg = get_wgs84_to_twd97_params()
% Parameterization:
cfg.a       = 6378137.0;                % Equatorial radius [unit in meter]
cfg.b       = 6356752.314245;           % Polar radius [unit in meter]
cfg.long0   = deg2rad(121);             % Central meridian of zone [unit in radian]
cfg.k0      = 0.9999;                   % Scale along long0
cfg.dx      = 250000;                   % Delta of shift in x-direction [unit in meter]
cfg.dy      = 0;                        % Delta of shift in y-direction [unit in meter]
end