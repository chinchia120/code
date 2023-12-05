clc;
clear all;
close all;

%%% ----- import file ----- %%%
filename = "group_data.xlsx";
[num, txt, raw] = xlsread(filename, 'A2: E23');

%%% ----- close net ----- %%%
net1 = ["B067", "B039", "B115", "BM1", "B067"];
net1_raw = cell([1, 5]);
cnt = 1;
for i = 1: 22
    for j = 1: 4
        if cell2mat(raw(i, 2)) == net1(j) & cell2mat(raw(i, 3)) == net1(j+1)
            net1_raw(cnt, :) = raw(i, :);
            cnt = cnt + 1;
        elseif cell2mat(raw(i, 2)) == net1(j+1) & cell2mat(raw(i, 3)) == net1(j)
            net1_raw(cnt, :) = raw(i, :);
            cnt = cnt + 1;
        end
    end 
end

%%% ----- chi-square test ----- %%%
std_before = 1.0;
X2_net1 = (7*(0.75037*10^(-3)*1000)^2)/(std_before^2)
X2_net2 = (7*(0.53731*10^(-3)*1000)^2)/(std_before^2)
X2_net3 = (11*(0.71303*10^(-3)*1000)^2)/(std_before^2)