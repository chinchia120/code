%% =============== Setup =============== %%
clc;
clear all;
close all;

%% =============== Q2 =============== %%
XA = 2017.358; % m
YA = 2538.658; % m
XB = 3009.888; % m
YB = 2101.241; % m

L_AC = 450; % m
HA_BAC = 230; % deg

phi_AB = rad2deg(atan((XB-XA) / (YB-YA))) + 180;
phi_AC = phi_AB + HA_BAC;
phi_CA = phi_AC + 180 - 360;

% XC = XA - L_AC * sind(360 - phi_AC);
% YC = YA + L_AC * cosd(360 - phi_AC);
XC = XA + L_AC * sind(phi_AC);
YC = YA + L_AC * cosd(phi_AC);

disp([sprintf('%%%% =============== Q2 =============== %%%%')]);
disp([sprintf('phi_CA = %06.3f(deg), XC = %7.3f(m), YC = %7.3f(m)\n', phi_CA, XC, YC)]);

%% =============== Q3 =============== %%
obs = [65.359, 65.352, 65.363, 65.354, 65.355, 65.352, 65.364, 65.358, 65.352]; % m

mean_obs = mean(obs);
scale = 99.975/100;
most_probable_value = mean_obs * scale;

obs = obs * scale;
v = 0;
for i = 1: size(obs, 2)
    v = v + (obs(i) - most_probable_value)^2;
end
std_obs = sqrt(v / (size(obs, 2)-1));
std_most_probable_value = sqrt(v / (size(obs, 2)-1) / (size(obs, 2)));

disp([sprintf('%%%% =============== Q3 =============== %%%%')]);
disp([sprintf('(a) Yes, the observation of 56.360(m) is blunder.')]);
disp([sprintf('(b) Most Probable Value = %5.3f(m)', most_probable_value)]);
disp([sprintf('(c) Presicion of Observation = %5.3f(m)', std_obs)]);
disp([sprintf('(d) Standard Deviation of Most Probable Value = %5.3f(m)\n', std_most_probable_value)]);

%% =============== Q4 =============== %%
syms L W H sigma_L sigma_W sigma_H;

V = L * W * H;

sigma_V = sqrt((diff(V, L, 1)^2)*sigma_L^2 + (diff(V, W, 1)^2)*sigma_W^2 + (diff(V, H, 1)^2)*sigma_H^2);
sigma_V = subs(sigma_V, [L, W, H], [50, 25, 5]);

sigma_V_sub_L = vpa(subs(sigma_V, [L, W, H, sigma_W, sigma_H], [50, 25, 5, 0.008, 0.004]));
sigma_V_sub_W = vpa(subs(sigma_V, [L, W, H, sigma_L, sigma_H], [50, 25, 5, 0.012, 0.004]));
sigma_V_sub_H = vpa(subs(sigma_V, [L, W, H, sigma_L, sigma_W], [50, 25, 5, 0.012, 0.008]));

disp([sprintf('%%%% =============== Q4 =============== %%%%')]);
disp([sprintf('If substitute L, sigma_v = %s(m)', sigma_V_sub_L)]);
disp([sprintf('If substitute W, sigma_v = %s(m)', sigma_V_sub_W)]);
disp([sprintf('If substitute H, sigma_v = %s(m)\n', sigma_V_sub_H)]);

%% =============== Q5 =============== %%
dH_AC = -4.026; % m //2.0km
dH_AD =  6.952; % m //2.0km
dH_DC = -3.042; % m //3.5km
dH_BC = -0.965; % m //2.0km
dH_BD =  2.031; % m //2.0km

dH_CA =  4.022; % m //2.5km
dH_DA = -6.948; % m //2.5km
dH_CD =  3.038; % m //4.0km
dH_CB =  0.961; % m //2.5km
dH_DB = -2.027; % m //2.5km

HA = 100;

close_BDC = dH_BD + dH_DC + dH_CB;
close_ADC = dH_AD + dH_DC + dH_CA;
dH_AC = -dH_AC;
dH_CA = -dH_CA;

% ===== AC ===== %
P_AC = 1/2.0;
P_CA = 1/2.5;
H_AC = (P_AC*dH_AC + P_CA*-dH_CA) / (P_AC+P_CA);

% ===== AD ===== %
P_AD = 1/2.0;
P_DA = 1/2.5;
H_AD = (P_AD*dH_AD + P_DA*-dH_DA) / (P_AD+P_DA);

% ===== BC ===== %
P_BC = 1/2.0;
P_CB = 1/2.5;
H_BC = (P_BC*dH_BC + P_CB*-dH_CB) / (P_BC+P_CB);

% ===== BD ===== %
P_BD = 1/2.0;
P_DB = 1/2.5;
H_BD = (P_BD*dH_BD + P_DB*-dH_DB) / (P_BD+P_DB);

% ===== DC ===== %
P_DC = 1/3.5;
P_CD = 1/4.0;
H_DC = (P_DC*dH_DC + P_CD*-dH_CD) / (P_DC+P_CD);

% ===== model ===== %
A = [ 0,  1,  0;
      0,  0,  1;
     -1,  1,  0;
     -1,  0,  1;
      0,  1, -1];

l = [H_AC + HA;
     H_AD + HA;
     H_BC     ;
     H_BD     ;
     H_DC     ;];

P = diag([(1/2+1/2.5), (1/2+1/2.5), (1/3.5+1/4), (1/2+1/2.5), (1/2+1/2.5)]);

dH = [H_AC;
      H_AD;
      H_BC;
      H_BD;
      H_DC;];

X_hat = inv(A'*P*A)*A'*P*l;
v_hat = A*X_hat - l;
dH_hat = dH + v_hat;
sigma_0_hat = sqrt((v_hat.'*P*v_hat)/(5-3));
sigma_HH_hat = sigma_0_hat^2 * (A.'*P*A);
sigma_HB = sqrt(sigma_HH_hat(1, 1));
sigma_HC = sqrt(sigma_HH_hat(2, 2));
sigma_HD = sqrt(sigma_HH_hat(3, 3));

disp([sprintf('%%%% =============== Q5 =============== %%%%')]);
disp([sprintf('(a) dH_AC and dH_CA')]);
disp([sprintf('(b) dH_AC = %4.3f(m), dH_AD = %4.3f(m), dH_BC = %4.3f(m), dH_BD = %4.3f(m), dH_DC = %4.3f(m)', H_AC, H_AD, H_BC, H_BD, H_DC)]);
disp([sprintf('(c) HB = %6.3f(m), HC = %6.3f(m), HD = %6.3f(m)', X_hat(1), X_hat(2), X_hat(3))]);
disp([sprintf('    v_H_AC = %4.3f(m), v_H_AD = %4.3f(m), v_H_BC = %4.3f(m), v_H_BD = %4.3f(m), v_H_DC = %4.3f(m)', v_hat(1), v_hat(2), v_hat(3), v_hat(4), v_hat(5))]);
disp([sprintf('    sigma_0_hat = %4.3f(m)', sigma_0_hat)]);
disp([sprintf('    dH_AC_hat = %4.3f(m), dH_AD_hat = %4.3f(m), dH_BC_hat = %4.3f(m), dH_BD_hat = %4.3f(m), dH_DC_hat = %4.3f(m)', dH_hat(1), dH_hat(2), dH_hat(3), dH_hat(4), dH_hat(5))]);
disp([sprintf('    sigma_HB = %4.3f(m), sigma_HC = %4.3f(m), sigma_HD = %4.3f(m)', sigma_HB, sigma_HC, sigma_HD)]);