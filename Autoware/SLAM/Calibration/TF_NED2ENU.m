%==========================================================================
% Function: NED to ENU Transform
% Developed & Coded by: Surachet Srinara, PhD candidate
% Date: April 28, 2021 - 12:06PM
% Copyright@Chet2021
%==========================================================================
function LLF_ENU_POSE_INTER = TF_NED2ENU(LLF_NED_POSE_INTER)
    [nRow, nCol] = size(LLF_NED_POSE_INTER);
    if (nCol == 19)
        LLF_ENU_POSE_INTER(1:nRow,1) = LLF_NED_POSE_INTER(:,1);                 % GPS time, t
        LLF_ENU_POSE_INTER(1:nRow,2) = LLF_NED_POSE_INTER(:,3);                 % East, x
        LLF_ENU_POSE_INTER(1:nRow,3) = LLF_NED_POSE_INTER(:,2);                 % North, y
        LLF_ENU_POSE_INTER(1:nRow,4) = -LLF_NED_POSE_INTER(:,4);                 % Height, z
        LLF_ENU_POSE_INTER(1:nRow,5) = LLF_NED_POSE_INTER(:,6);                 % V_E, x
        LLF_ENU_POSE_INTER(1:nRow,6) = LLF_NED_POSE_INTER(:,5);                 % V_N, y
        LLF_ENU_POSE_INTER(1:nRow,7) = -LLF_NED_POSE_INTER(:,7);                % V_U, z
        LLF_ENU_POSE_INTER(1:nRow,8) = LLF_NED_POSE_INTER(:,9);                 % Pitch, x
        LLF_ENU_POSE_INTER(1:nRow,9) = LLF_NED_POSE_INTER(:,8);                 % Roll, y
        LLF_ENU_POSE_INTER(1:nRow,10) = -LLF_NED_POSE_INTER(:,10);              % Heading, z
        LLF_ENU_POSE_INTER(1:nRow,11) = LLF_NED_POSE_INTER(:,12);               % E-position STD, x
        LLF_ENU_POSE_INTER(1:nRow,12) = LLF_NED_POSE_INTER(:,11);               % N-position STD, y
        LLF_ENU_POSE_INTER(1:nRow,13) = LLF_NED_POSE_INTER(:,13);               % U-position STD, z
        LLF_ENU_POSE_INTER(1:nRow,14) = LLF_NED_POSE_INTER(:,15);               % E-velocity STD, x
        LLF_ENU_POSE_INTER(1:nRow,15) = LLF_NED_POSE_INTER(:,14);               % N-velocity STD, y
        LLF_ENU_POSE_INTER(1:nRow,16) = LLF_NED_POSE_INTER(:,16);               % U-velocity STD, z
        LLF_ENU_POSE_INTER(1:nRow,17) = LLF_NED_POSE_INTER(:,18);               % Pitch angle STD, x
        LLF_ENU_POSE_INTER(1:nRow,18) = LLF_NED_POSE_INTER(:,17);               % Roll angle STD, y
        LLF_ENU_POSE_INTER(1:nRow,19) = LLF_NED_POSE_INTER(:,19);               % Heading angle STD, z
        
    elseif (nCol == 10)
        LLF_ENU_POSE_INTER(1:nRow,1) = LLF_NED_POSE_INTER(:,1);                 % GPS time, t
        LLF_ENU_POSE_INTER(1:nRow,2) = LLF_NED_POSE_INTER(:,3);                 % East, x
        LLF_ENU_POSE_INTER(1:nRow,3) = LLF_NED_POSE_INTER(:,2);                 % North, y
        LLF_ENU_POSE_INTER(1:nRow,4) = -LLF_NED_POSE_INTER(:,4);                 % Height, z
        LLF_ENU_POSE_INTER(1:nRow,5) = LLF_NED_POSE_INTER(:,6);                 % V_E, x
        LLF_ENU_POSE_INTER(1:nRow,6) = LLF_NED_POSE_INTER(:,5);                 % V_N, y
        LLF_ENU_POSE_INTER(1:nRow,7) = LLF_NED_POSE_INTER(:,7);                 % V_U, z
        LLF_ENU_POSE_INTER(1:nRow,8) = LLF_NED_POSE_INTER(:,9);                 % Pitch, x
        LLF_ENU_POSE_INTER(1:nRow,9) = LLF_NED_POSE_INTER(:,8);                 % Roll, y
        LLF_ENU_POSE_INTER(1:nRow,10) = -LLF_NED_POSE_INTER(:,10);              % Heading, z
        
    elseif (nCol == 13) % GNSS solution
        LLF_ENU_POSE_INTER(1:nRow,1) = LLF_NED_POSE_INTER(:,1);                 % GPS time, t
        LLF_ENU_POSE_INTER(1:nRow,2) = LLF_NED_POSE_INTER(:,3);                 % East, x
        LLF_ENU_POSE_INTER(1:nRow,3) = LLF_NED_POSE_INTER(:,2);                 % North, y
        LLF_ENU_POSE_INTER(1:nRow,4) = -LLF_NED_POSE_INTER(:,4);                 % Height, z
        LLF_ENU_POSE_INTER(1:nRow,5) = LLF_NED_POSE_INTER(:,6);                 % V_E, x
        LLF_ENU_POSE_INTER(1:nRow,6) = LLF_NED_POSE_INTER(:,5);                 % V_N, y
        LLF_ENU_POSE_INTER(1:nRow,7) = LLF_NED_POSE_INTER(:,7);                 % V_U, z
        LLF_ENU_POSE_INTER(1:nRow,8) = LLF_NED_POSE_INTER(:,9);                 % E-position STD, x
        LLF_ENU_POSE_INTER(1:nRow,9) = LLF_NED_POSE_INTER(:,8);                 % N-position STD, y
        LLF_ENU_POSE_INTER(1:nRow,10) = LLF_NED_POSE_INTER(:,10);               % U-position STD, z
        LLF_ENU_POSE_INTER(1:nRow,11) = LLF_NED_POSE_INTER(:,12);               % E-velocity STD, x
        LLF_ENU_POSE_INTER(1:nRow,12) = LLF_NED_POSE_INTER(:,11);               % N-velocity STD, y
        LLF_ENU_POSE_INTER(1:nRow,13) = LLF_NED_POSE_INTER(:,13);               % U-velocity STD, z

    end
end