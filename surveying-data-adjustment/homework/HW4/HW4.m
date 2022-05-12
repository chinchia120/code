clc;
clear all;
format long

%Q3
L=[dms2deg([70 00 00]) ; dms2deg([50 00 00]) ; dms2deg([-119 -59 -30])];
L_dms=rad2dms(deg2rad(L));
A=[1 0 ; 0 1 ; -1 -1];
AT=transpose(A);
X=inv(AT*A)*(AT*L);
X_dms=rad2dms(deg2rad(X));
% V=A*X-L;
% V_dms=rad2dms(deg2rad(V));
% V_len=[V_dms(1,3) ; V_dms(2,3) ; V_dms(3,3)]/206265;
% tmp=transpose(V_len)*V_len*inv(AT*A);
tmp=[(20/206265)^2 0 ; 0 (20/206265)^2];
AC=dms2deg([336 51 52]);
J=[500*(cosd(X(1,1))/sind(X(2,1))*sind(AC)+sind(X(1,1))/sind(X(2,1))*cosd(AC)) 500*(-sind(X(1,1))*cosd(X(2,1))/sind(X(2,1))/sind(X(2,1))*sind(AC)+sind(X(1,1))/sind(X(2,1))*cosd(AC)) ;
   500*(cosd(X(1,1))/sind(X(2,1))*cosd(AC)-sind(X(1,1))/sind(X(2,1))*sind(AC)) 500*(-sind(X(1,1))*cosd(X(2,1))/sind(X(2,1))/sind(X(2,1))*cosd(AC)-sind(X(1,1))/sind(X(2,1))*sind(AC))];
tmp2=J*[tmp(1,1) 0 ; 0 tmp(2,2)]*transpose(J);
SD=sqrt(tmp2);

%Q5
length=[50.2 ; 120.5 ; 200.0];
angle=[dms2rad([00 27 00]) ; dms2rad([01 05 00]) ; dms2rad([01 47 00])];
L=length./angle;
A=[1 ; 1 ; 1];
AT=transpose(A);
P=diag([1/((1/angle(1,1)*0.4)^2+(length(1,1)/angle(1,1)^2*dms2rad([00 0.2 00]))^2) 1/((1/angle(2,1)*0.4)^2+(length(2,1)/angle(2,1)^2*dms2rad([00 0.2 00]))^2) 1/((1/angle(3,1)*0.4)^2+(length(3,1)/angle(3,1)^2*dms2rad([00 0.2 00]))^2)]);
X=inv(AT*P*A)*(AT*P*L);
V=A*X-L;
tmp=transpose(V)*P*V/(3-1)*inv(AT*P*A);
SD=sqrt(tmp);

%Q10
syms S1;
syms S2;
syms S3;
syms S4;
syms S5;
syms S6;
syms S7;
syms S8;
syms S9;
syms S10;

syms SD1;
syms SD2;
syms SD3;
syms SD4;
syms SD5;
syms SD6;
syms SD7;
syms SD8;
syms SD9;
syms SD10;

L=[S1 ; S2 ; S3 ; S4 ; S5 ; S6 ; S7 ; S8 ; S9 ; S10];
P=diag([SD1 SD2 SD3 SD4 SD5 SD6 SD7 SD8 SD9 SD10]);
A=[-1 1 0 0 0 ; 0 -1 1 0 0 ; 0 0 -1 1 0 ; 0 0 0 -1 1 ; -1 0 1 0 0 ; 0 0 -1 0 1 ; -1 0 0 1 0 ; 0 -1 0 1 0 ; 0 -1 0 0 1 ; -1 0 0 0 1];
AT=transpose(A);
X=inv(AT*P*A)*(AT*P*L);

%Q18
L=[411.73 ; -694.61 ; 816.38 ; 1671.17 ; -22.27 ; 82.13 ; 214.75 ; -852.25];
A=[350 -238 1 0 ; -750 262 1 0 ; 850 -138 1 0 ; 1350 -1538 1 0 ; -238 -350 0 1 ; 262 750 0 1 ; -138 -850 0 1 ; -1538 -1350 0 1];
AT=transpose(A);
X=inv(AT*A)*(AT*L);
V=A*X-L;
SD_0=sqrt((transpose(V)*V)/4);
matrix=SD_0^2*inv(AT*A);
sd=sqrt(matrix);
coff_ab=matrix(1,2)/(sd(1,1)*sd(2,2));
coff_ae=matrix(1,3)/(sd(1,1)*sd(3,3));
coff_af=matrix(1,4)/(sd(1,1)*sd(4,4));
coff_be=matrix(2,3)/(sd(2,2)*sd(3,3));
coff_bf=matrix(2,4)/(sd(2,2)*sd(4,4));
coff_ef=matrix(3,4)/(sd(3,3)*sd(4,4));

%Q22
deg1=0;
deg2=0;
deg3=0;
deg4=40;
deg5=10;
deg6=25;
deg7=-40;
L=[deg1 ; deg2 ; deg3 ; deg4 ; deg5 ; deg6 ; deg7];
A=[1 0 0 ; 0 1 0 ; 0 0 1 ; 1 1 0 ; 0 1 1 ; 1 1 1 ; -1 -1 -1];
AT=transpose(A);
X=inv(AT*A)*(AT*L);
V=A*X-L;
SD=sqrt(transpose(V)*V/4*A*inv(AT*A)*AT);

%Q23
L=[10.008+50 ; 19.987 ; 10.012 ; -15.004 ; -4.992 ; 10.010 ; 4.998-70];
A=[0 0 0 1 ; 1 0 0 -1 ; 0 0 1 -1 ; -1 1 0 0 ; 0 1 -1 0 ; 1 0 -1 0 ; 0 -1 0 0];
AT=transpose(A);
X=inv(AT*A)*(AT*L);
V=A*X-L;
tmp=transpose(V)*V/3*inv(AT*A);
SD=sqrt(tmp);
SD_PR=sqrt([-1 1]*([tmp(1,1) tmp(1,3) ; tmp(1,3) tmp(3,3)])*[-1 ; 1]);
SD_QS=sqrt([-1 1]*([tmp(2,2) tmp(2,4) ; tmp(2,4) tmp(4,4)])*[-1 ; 1]);

%Q24
L=[72.840+110.394 ; -16.772+200.052 ; 14.802+168.478 ; 39.774+143.506];
A=[1 ; 1 ; 1 ; 1];
AT=transpose(A);
P=diag([1/25 ; 1/18 ; 1/15 ; 1/12]);
X=inv(AT*P*A)*(AT*P*L);
V=A*X-L;
sigma0=sqrt(transpose(V)*P*V/3);
tmp=transpose(V)*P*V/3*inv(AT*P*A);
SD=sqrt(tmp);