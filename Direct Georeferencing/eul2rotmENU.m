%% =============== Function of eul2rotmENU =============== %%
function rot = eul2rotmENU(eulAngIn)
      % == Euler's angles around the XYZ (Right-Fwd-Up) axes:
      omeX = eulAngIn(1);                   % Pitch angle (X-RIGHT)
      phiY = eulAngIn(2);                   % Roll angle (Y-FWD)
      kapZ = eulAngIn(3);                   % Heading angle (Z-Up)

      % == Rotation matrix of X-Right axis:
      Rx = [1, 0, 0; 
            0,  cosd(omeX), sind(omeX);
            0, -sind(omeX), cosd(omeX)];

      % == Rotation matrix of Y-Fwd axis:
      Ry = [cosd(phiY), 0, -sind(phiY); 
            0, 1, 0;
            sind(phiY), 0, cosd(phiY)];  

      % == Rotation matrix of Z-Up axis:
      Rz = [ cosd(kapZ), sind(kapZ), 0; 
            -sind(kapZ), cosd(kapZ), 0;
            0, 0, 1];    

      % == Output:
      rot.matrix = (Rx*(Ry*Rz));    
      eulAngOut = getEulerFromRotm(rot.matrix, eulAngIn);
      rot.euler(1: 3, 1) = eulAngOut;
end