%%% Q1 %%%

%% Initial Data %%%
syms a b c d e f g h k X Y sigma_x sigma_y;
U = a*X + b*Y + c;
V = d*X + e*Y + f;
Z = g*U + h*V + k;

sigma_XX = [sigma_x^2, 0; 0, sigma_y^2];

%% (a) %%
J = [diff(Z, X), diff(Z, Y)];
sigma_z = sqrt(J*sigma_XX*J.');

%% (b) %%
J = [diff(U, X), diff(U, Y); diff(V, X), diff(V, Y); diff(Z, X), diff(Z, Y)];
cov_UVZ = J*sigma_XX*J.';

%% (c) %%
J = [diff(Z, X), diff(Z, Y); diff(X, X), diff(X, Y)];
cov_ZX = J*sigma_XX*J.';

J = [diff(Z, X), diff(Z, Y); diff(Y, X), diff(Y, Y)];
cov_ZY = J*sigma_XX*J.';

%%% Q2 %%%

%% (a) %%
A1 = 300.25; 
A2 = 300.35; 

avg_A = (A1+A2)/2;
std_A = sqrt(((A1-avg_A)^2+(A2-avg_A)^2)/(2-1));

%% (b) %%
SEM = std_A/sqrt(2);

%%% Q4 %%% 

%% (a) %%
sigma = 7*sqrt(2);

%% (b) %%
sigma_fb = sqrt(2*sigma^2);

%% (c) %%
sigma_h = sqrt(2*(1/2)^2*sigma^2);

%%% Q5 %%%

%% (b) %%
sigma_b = 3 * sqrt(4/5);
sigma_a = 2 * sigma_b;

%% (c) %%
sigma_A = sigma_a / sqrt(6);
sigma_B = sigma_b / sqrt(4);