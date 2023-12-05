clc;
clear;
format long;

%Q1
B=[1 1 1];
W=[-20];
P=diag([1 2 4]);
V=inv(P)*transpose(B)*inv(B*inv(P)*transpose(B))*W;
%V_dms = rad2dms(deg2rad(V))
SD_0 = sqrt(V.'*P*V/1);