%% --------------- Setup --------------- %%
clc;
clear all;
close all;

%% --------------- PVAF --------------- %%
r = 10; % percent of year
n = 6; % number of year

PVAF = (1-(1/(1+r*0.01)^n)) / (r*0.01);
disp([sprintf('PVAF = %.3f\n', PVAF)]);

%% --------------- FVAF --------------- %%
r = 2.2/12; % percent of month
n = 80*12; % number of month

FVAF = ((1+r*0.01)^n-1)/(r*0.01);
disp([sprintf('FVAF = %.3f\n', FVAF)]);

%% --------------- Valus of Bond --------------- %%
Coupon_Rate = 0; % percent of year
Face_Value = 1000; % dolar
Yield_to_Maturity = 8.8; % percent of year
n = 10; % number of year
Coupon_Payment = Face_Value * (Coupon_Rate*0.01) * 1/2; % dolar / half of year

PVAF = (1-(1/(1+Yield_to_Maturity/2*0.01)^(n*2))) / (Yield_to_Maturity/2*0.01);
VB = Coupon_Payment*PVAF + Face_Value/(1+Yield_to_Maturity/2*0.01)^(n*2);
disp([, sprintf('VB = %.4f\n', VB)]);

%% --------------- Holding Period Return --------------- %%
Coupon_Rate = 8; % percent of year
Face_Value = 1000; % dolar
Yield_to_Maturity = 10; % percent of year
n1 = 20; % number of year
n2 = 19; % number of year
Coupon_Payment = Face_Value * (Coupon_Rate*0.01) * 1/2; % dolar / half of year

PVAF = (1-(1/(1+Yield_to_Maturity/2*0.01)^(n1*2))) / (Yield_to_Maturity/2*0.01);
VB1 = Coupon_Payment*PVAF + Face_Value/(1+Yield_to_Maturity/2*0.01)^(n1*2);

PVAF = (1-(1/(1+Yield_to_Maturity/2*0.01)^(n2*2))) / (Yield_to_Maturity/2*0.01);
VB2 = Coupon_Payment*PVAF + Face_Value/(1+Yield_to_Maturity/2*0.01)^(n2*2);

HPR = (VB2 - VB1)/VB1 + Coupon_Payment*2/VB1; 
disp([sprintf('VB1 = %.4f, VB2 = %.4f, HPR = %2.2f%%\n', VB1, VB2, HPR*100)]);

%% --------------- Current Yield --------------- %%
Coupon_Rate = 8; % percent of year
Face_Value = 1000; % dolar
Yield_to_Maturity = 6; % percent of year
n = 20; % number of year
Coupon_Payment = Face_Value * (Coupon_Rate*0.01) * 1/2; % dolar / half of year

PVAF = (1-(1/(1+Yield_to_Maturity/2*0.01)^(n*2))) / (Yield_to_Maturity/2*0.01);
VB = Coupon_Payment*PVAF + Face_Value/(1+Yield_to_Maturity/2*0.01)^(n*2);
CY = Coupon_Payment*2 / VB;
disp([sprintf('VB = %.4f, CY = %2.2f%%\n', VB, CY*100)]);
