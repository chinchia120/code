% Get ENU-based Navigation Solution from Smoothed Navigation Data
% Developed & Coded by: Surachet Srinara
% June 4, 2022

function nav = getNavSolENU(f_nav, rowNum)
% == Get Euler angles from an input nav. file:
omegaX          = (f_nav(rowNum, 8));     % Pitch angle (X-RIGHT)
phiY            = (f_nav(rowNum, 9));     % Roll angle (Y-FWD)
kappaZ          = (f_nav(rowNum, 10));    % -Heading angle (Z-UP)
eulerAng        = [omegaX; phiY; kappaZ]; % [roll,pitch,heading]

% == Rotations follow the right-handed rule sequencing with ZYX:
% rot             = eul2rotmENU(eulerAng); % C_nb
rot = eul2Cnb([phiY,omegaX,kappaZ]); % C_nb (n-frame is ENU to RFU)
nav.rotm        = rot.matrix;

% == Translation vector:
nav.tran        = f_nav(rowNum, 2:4);

% == Velocity:
nav.vel         = f_nav(rowNum, 5:7);

% == Attitude:
nav.att         = rot.euler;
% if norm(rot.euler-eulerAng) > (1/100000)
%     nav.att = eulerAng;
% end

end