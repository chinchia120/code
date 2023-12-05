%% =============== Function of Cbn2eul =============== %%
function PSI_nb = Cbn2eul(Cbn)
    % body frame: R-F-U
    % n-frame: E-N-U
    Cnb = Cbn'; % C_nb

    pitch = asind(Cnb(2, 3));
    heading = atan2d(-Cnb(2, 1),Cnb(2, 2));
    roll = atan2d(-Cnb(1, 3),Cnb(3, 3));

    PSI_nb = [roll, pitch, -heading];
end
