%==========================================================================
% Function: INS/GNSS-based POS Interpolation
% Developed & Coded by: Surachet Srinara, PhD candidate
% Date: July 8, 2021 - 11:50AM
% Copyright@Chet2021
%==========================================================================
function inter_POS = interpolatePOS(lidarTime, f_llf_ins, ins_row)
[nRow, nCol] = size(f_llf_ins);
    
% INS-Time synchronization and interpolation:    
t_ins0 = f_llf_ins(ins_row-1,1);
t_ins1 = f_llf_ins(ins_row,1);
inter_POS = zeros(1, nCol);

for n = 1:nCol
    if (n == 1)
        inter_POS(1,n) = lidarTime;              % GPS time stamp
    else
        inter_POS(1,n) = f_llf_ins(ins_row,n) - ...
            ((f_llf_ins(ins_row,n) - f_llf_ins(ins_row-1,n))/((t_ins1 - t_ins0)/(t_ins1 - lidarTime)));
    end
end
