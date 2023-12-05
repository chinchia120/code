clc;
clear all;
close all;
%% ---------- output file ---------- %%
file = fopen('HW4_output.txt', 'w');
fprintf(file, '%%%% ---------- Q1 ---------- %%%%\n');

%% ---------- Q1 ---------- %%
v = [7; -4; 2; -18; -7];
r = [0.36, 0.25, 0.03, 0.48, 0.64];
std = 5; 
z = 1.96;

for i = 1: size(v, 1)
    w(i, 1) = {abs(v(i))/(std*sqrt(r(i)))};
    
    if(cell2mat(w(i, 1)) > z)
        str = "fail";
    else
        str = "pass";
    end
    w(i, 2) = {str};
    fprintf(file, 'Case %d: v = %3d, w = %3.2f, check = %s\n', i, v(i), cell2mat(w(i, 1)), cell2sym(w(i, 2)));
end

%% ---------- Q2 ---------- %%
QP = [0.71, 0.59, 1, 0.59, 0.71, 0.35, 0.35, 0.71, 0];
sum_QP_diag = sum(QP);
fprintf(file, '\n%%%% ---------- Q2 ---------- %%%%\n');
fprintf(file, '(a): 1. 是一種冪等矩陣\n\t 2. 是一種秩缺方陣\n\t 3. 其跡等於平差系統的總多於觀測數 \n\t ');
fprintf(file, '4. 若為不相關觀測，即P矩陣為對角矩陣，則QP矩陣的第i個主對角線元素稱為第i個觀測量的多餘分量，或稱\n\t ');
fprintf(file, '   為局部多於觀測數記為ri，且0 <= ri <= 1。在這種情況下，當觀測量的多餘分量ri = 0時，表示該觀測\n\t ');
fprintf(file, '   量可能為必要觀測量，若配合其所屬行列的QP個元素值均為0，則必為必要觀測量；若ri = 1表示該觀測量\n\t ');
fprintf(file, '   完全是多餘的。QP矩陣和設計矩陣A及觀測量的權矩陣P有關，故ri和觀測精度P及幾何結構A有關\n');
fprintf(file, '(b): 1. i觀測量若有粗差，會完全無法被偵測出來\n\t ');
fprintf(file, '2. QP(9, :) = 0\n');
fprintf(file, '(c): 1. c觀測量若有粗差，會完全反應在該觀測量本身的改正數上\n\t ');
fprintf(file, '2. QP(3, 3) = 1 && QP(3, !3:) = 0\n');
fprintf(file, '(d): 1. c和i觀測量若有粗差，對其他觀測量的改正數沒有影響\n\t ');
fprintf(file, '2. QP(3, !3:) = 0, QP(9, !9:) = 0\n');
fprintf(file, '(e): 1. i觀測量為必要觀測量\n\t ');
fprintf(file, '2. QP(9, :) = 0\n');
fprintf(file, '(f): 1. 此水準網平差系統有4個未知參數個數\n\t ');
fprintf(file, '2. sum_QP_diag = %3.2f, %d - %3.2f = %3.2f, round(%3.2f) = %3.2f \n', sum_QP_diag, size(QP, 2), sum_QP_diag, size(QP, 2)-sum_QP_diag, size(QP, 2)-sum_QP_diag, round(size(QP, 2)-sum_QP_diag));

