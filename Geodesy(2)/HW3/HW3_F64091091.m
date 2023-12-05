clc;
clear;
close;
file = fopen('HW3_F64091091_output.txt', 'w');
%% ---------- Q1 ---------- %%
W0 = 62636850.0;
WA = 62636860.0;
CA = W0 - WA;

lat = 0.0;
lon = 0.0;
N = 10.0;
g = 9.780;

a = 6378137.0;
b = 6356752.314;
f = 0.00335281068118;
m = 0.00344978600308;
e1 = sqrt(a^2-b^2) / a;

% ----- Dynamic Height ----- %
re = 9.7803267715;
k = 0.001931851353;
r0_45 = re*(1+k*sind(45)^2)/sqrt(1-e1^2*sind(45)^2);
dynamic_height = CA / r0_45;

% ----- Orthometric Height ----- %
orthometric_height = (-g+sqrt(g^2+4*0.0424*CA))/(2*0.0424);

% ----- Normal Height ----- %
r_bar = r0_45 - (0.15438-0.00022*sind(lat)^2)*orthometric_height + 1.8*10^(-7)*orthometric_height^2*10^(-5);
normal_height = CA / r_bar;

% ----- Ellipsoidal Height ----- %
ellipsoidal_height = orthometric_height + N;

%% ---------- Q7 ---------- %%
Sij = 1000;
Zij = 88; 
ii = 1.5;
ti = 1.2;
kij = 1.3;
Dij = Sij*sind(Zij);
deflection_of_the_vertical = 1*(60/206265);

lat = 30;
azimuth = 30;
M = (a*(1-(e1^2)))/((1-(e1^2)*(sin(lat)^2))^(3/2));
N = a/((1-(e1^2)*(sin(lat)^2))^(1/2));
R = 1/(cosd(azimuth)^2/M + sind(azimuth)^2/N);

hij = Sij*cosd(Zij) + (Dij^2)/(2*R) - (kij*Dij^2)/(2*R) - deflection_of_the_vertical*Dij + ii - ti;

%% ---------- Q9 ---------- %%
R = 6378137;
a = 6.7*10^(-8);
c = 5*(60/206265);
dT = 20.0-20.1;
foresight_dis = 30.00;
backsight_dis = 30.10;
foresight_red = 1.200;
backsight_red = 1.100;
dn = backsight_red-foresight_red;
S = (foresight_dis+backsight_dis)/2;

% ----- Correction of Refraction Errors ----- %
RE = -a*(S^2)*dT*dn;
dn_refraction = dn + RE;

% ------ Correction of Earth's Curvature ----- %
db = (backsight_dis^2)/(2*R);
df = (foresight_dis^2)/(2*R);
dn_curvature = dn - (db-df);

% ------ Correction of Collimation Errors ----- %
dn_collimation = backsight_red - foresight_red - c*(backsight_dis-foresight_dis);

%% ---------- Output ---------- %%
fprintf(file, '%%%% ---------- Q1 ---------- %%%%\n');
fprintf(file, 'Dynamic Height = %5.5f (m)\n', dynamic_height);
fprintf(file, 'Orthometric Height = %5.5f (m)\n', orthometric_height);
fprintf(file, 'Normal Height = %5.5f (m)\n', normal_height);
fprintf(file, 'Ellipsoidal Height = %5.5f (m)\n\n', ellipsoidal_height);

fprintf(file, '%%%% ---------- Q7 ---------- %%%%\n');
fprintf(file, 'Dij = %4.4f (m)\n', Dij);
fprintf(file, 'M = %4.4f (m)\n', M);
fprintf(file, 'N = %4.4f (m)\n', N);
fprintf(file, 'R = %4.4f (m)\n', R);
fprintf(file, 'hij = %4.4f (m)\n\n', hij);

fprintf(file, '%%%% ---------- Q9 ---------- %%%%\n');
fprintf(file, 'dn = %10.10f (m)\n\n', dn);

fprintf(file, '%% --- Correction of Refraction Errors --- %%\n');
fprintf(file, 'RE = %4.4e (m)\n', RE);
fprintf(file, 'dn_refraction = %10.10f (m)\n\n', dn_refraction);

fprintf(file, '%% --- Correction of Earth''s Curvature --- %%\n');
fprintf(file, 'db = %4.4e (m)\n', db);
fprintf(file, 'df = %4.4e (m)\n', df);
fprintf(file, 'dn_curvature = %10.10f (m)\n\n', dn_curvature);

fprintf(file, '%% --- Correction of Collimation Errors --- %%\n');
fprintf(file, 'dn_collimation = %10.10f (m)\n\n', dn_collimation);