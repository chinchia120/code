%% =============== Function of eul2Cnb =============== %%
function rot = eul2Cnb(PSI_nb)
    % Rotate Z
    C3 = [cosd(PSI_nb(3)), sind(PSI_nb(3)), 0;
         -sind(PSI_nb(3)), cosd(PSI_nb(3)), 0;
          0, 0, 1];
    % Rotate X
    C2 = [1, 0, 0;
          0,  cosd(PSI_nb(2)), sind(PSI_nb(2));
          0, -sind(PSI_nb(2)), cosd(PSI_nb(2))];
    % Rotate Y
    C1 = [cosd(PSI_nb(1)), 0, -sind(PSI_nb(1));
          0, 1, 0;
          sind(PSI_nb(1)), 0,  cosd(PSI_nb(1))];

    rot.matrix = (C1*C2*C3);
    rot.euler = Cbn2eul(rot.matrix');  %[roll, pitch, heading]
end
