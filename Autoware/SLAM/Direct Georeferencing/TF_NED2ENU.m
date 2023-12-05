%% =============== Function of TF_NED2ENU =============== %%
function LLF_ENU_POSE_INTER = TF_NED2ENU(LLF_NED_POSE_INTER)
    [nRow, nCol] = size(LLF_NED_POSE_INTER);
    if (nCol == 19)
        LLF_ENU_POSE_INTER(1: nRow, 01) =  LLF_NED_POSE_INTER(:, 01);         % GPS time, t
        LLF_ENU_POSE_INTER(1: nRow, 02) =  LLF_NED_POSE_INTER(:, 03);         % East, x
        LLF_ENU_POSE_INTER(1: nRow, 03) =  LLF_NED_POSE_INTER(:, 02);         % North, y
        LLF_ENU_POSE_INTER(1: nRow, 04) = -LLF_NED_POSE_INTER(:, 04);         % Height, z
        LLF_ENU_POSE_INTER(1: nRow, 05) =  LLF_NED_POSE_INTER(:, 06);         % V_E, x
        LLF_ENU_POSE_INTER(1: nRow, 06) =  LLF_NED_POSE_INTER(:, 05);         % V_N, y
        LLF_ENU_POSE_INTER(1: nRow, 07) = -LLF_NED_POSE_INTER(:, 07);         % V_U, z
        LLF_ENU_POSE_INTER(1: nRow, 08) =  LLF_NED_POSE_INTER(:, 09);         % Pitch, x
        LLF_ENU_POSE_INTER(1: nRow, 09) =  LLF_NED_POSE_INTER(:, 08);         % Roll, y
        LLF_ENU_POSE_INTER(1: nRow, 10) = -LLF_NED_POSE_INTER(:, 10);         % Heading, z
        LLF_ENU_POSE_INTER(1: nRow, 11) =  LLF_NED_POSE_INTER(:, 12);         % E-position STD, x
        LLF_ENU_POSE_INTER(1: nRow, 12) =  LLF_NED_POSE_INTER(:, 11);         % N-position STD, y
        LLF_ENU_POSE_INTER(1: nRow, 13) =  LLF_NED_POSE_INTER(:, 13);         % U-position STD, z
        LLF_ENU_POSE_INTER(1: nRow, 14) =  LLF_NED_POSE_INTER(:, 15);         % E-velocity STD, x
        LLF_ENU_POSE_INTER(1: nRow, 15) =  LLF_NED_POSE_INTER(:, 14);         % N-velocity STD, y
        LLF_ENU_POSE_INTER(1: nRow, 16) =  LLF_NED_POSE_INTER(:, 16);         % U-velocity STD, z
        LLF_ENU_POSE_INTER(1: nRow, 17) =  LLF_NED_POSE_INTER(:, 18);         % Pitch angle STD, x
        LLF_ENU_POSE_INTER(1: nRow, 18) =  LLF_NED_POSE_INTER(:, 17);         % Roll angle STD, y
        LLF_ENU_POSE_INTER(1: nRow, 19) =  LLF_NED_POSE_INTER(:, 19);         % Heading angle STD, z
        
    elseif (nCol == 10)
        LLF_ENU_POSE_INTER(1: nRow, 01) =  LLF_NED_POSE_INTER(:, 01);         % GPS time, t
        LLF_ENU_POSE_INTER(1: nRow, 02) =  LLF_NED_POSE_INTER(:, 03);         % East, x
        LLF_ENU_POSE_INTER(1: nRow, 03) =  LLF_NED_POSE_INTER(:, 02);         % North, y
        LLF_ENU_POSE_INTER(1: nRow, 04) = -LLF_NED_POSE_INTER(:, 04);         % Height, z
        LLF_ENU_POSE_INTER(1: nRow, 05) =  LLF_NED_POSE_INTER(:, 06);         % V_E, x
        LLF_ENU_POSE_INTER(1: nRow, 06) =  LLF_NED_POSE_INTER(:, 05);         % V_N, y
        LLF_ENU_POSE_INTER(1: nRow, 07) =  LLF_NED_POSE_INTER(:, 07);         % V_U, z
        LLF_ENU_POSE_INTER(1: nRow, 08) =  LLF_NED_POSE_INTER(:, 09);         % Pitch, x
        LLF_ENU_POSE_INTER(1: nRow, 09) =  LLF_NED_POSE_INTER(:, 08);         % Roll, y
        LLF_ENU_POSE_INTER(1: nRow, 10) = -LLF_NED_POSE_INTER(:, 10);         % Heading, z
        
    elseif (nCol == 13) % GNSS solution
        LLF_ENU_POSE_INTER(1: nRow, 01) =  LLF_NED_POSE_INTER(:, 01);         % GPS time, t
        LLF_ENU_POSE_INTER(1: nRow, 02) =  LLF_NED_POSE_INTER(:, 03);         % East, x
        LLF_ENU_POSE_INTER(1: nRow, 03) =  LLF_NED_POSE_INTER(:, 02);         % North, y
        LLF_ENU_POSE_INTER(1: nRow, 04) = -LLF_NED_POSE_INTER(:, 04);         % Height, z
        LLF_ENU_POSE_INTER(1: nRow, 05) =  LLF_NED_POSE_INTER(:, 06);         % V_E, x
        LLF_ENU_POSE_INTER(1: nRow, 06) =  LLF_NED_POSE_INTER(:, 05);         % V_N, y
        LLF_ENU_POSE_INTER(1: nRow, 07) =  LLF_NED_POSE_INTER(:, 07);         % V_U, z
        LLF_ENU_POSE_INTER(1: nRow, 08) =  LLF_NED_POSE_INTER(:, 09);         % E-position STD, x
        LLF_ENU_POSE_INTER(1: nRow, 09) =  LLF_NED_POSE_INTER(:, 08);         % N-position STD, y
        LLF_ENU_POSE_INTER(1: nRow, 10) =  LLF_NED_POSE_INTER(:, 10);         % U-position STD, z
        LLF_ENU_POSE_INTER(1: nRow, 11) =  LLF_NED_POSE_INTER(:, 12);         % E-velocity STD, x
        LLF_ENU_POSE_INTER(1: nRow, 12) =  LLF_NED_POSE_INTER(:, 11);         % N-velocity STD, y
        LLF_ENU_POSE_INTER(1: nRow, 13) =  LLF_NED_POSE_INTER(:, 13);         % U-velocity STD, z
    end
end