% Get Euler angles from rotation matrix
% Developed by Surachet Srinara
% June 5, 2022

function eulAngOut = getEulerFromRotm(Rotm, eulAngIn)
% Get Euler angles from an input of rotation matrix:
eulAngOut           = rotm2eul(Rotm', 'ZYX')*180/pi;     
eulAngOut           = [eulAngOut(3); eulAngOut(2); eulAngOut(1)];

% Check an output:
% if (eulAngIn(3) < 0) && (eulAngOut(3) > 0) && (abs(eulAngIn(3)) > 90) && (abs(eulAngOut(3)) > 90)
%     preHeading      = eulAngOut(3) - 360;
%     eulAngOut(3,1)  = preHeading;
% end
% if (eulAngIn(3) > 0) && (eulAngOut(3) < 0) && (abs(eulAngIn(3)) > 90) && (abs(eulAngOut(3)) > 90)
%     preHeading      = 360 + eulAngOut(3);
%     eulAngOut(3,1)  = preHeading;
% end
% if norm(abs(eulAngOut(3))-abs(eulAngIn(3))) >= 5 
%     preHeading      = eulAngOut(3) - 360;
%     eulAngOut(3,1)  = preHeading;
% end
headingOffset = abs(eulAngOut(3)) - abs(eulAngIn(3));
end