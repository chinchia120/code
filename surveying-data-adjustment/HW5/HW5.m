clc;
clear all;
format long;

%Q1
B=[1 1 1];
W=[-18];
P=diag([1/2 1/3 1/4]);
V=inv(P)*transpose(B)*inv(B*inv(P)*transpose(B))*W;

%Q3
B=[1 1 1];
W=[-30];
P=diag([1 1 1]);
L=[dms2deg([70 00 00]) ; dms2deg([50 00 00]) ; dms2deg([60 00 30])];
V=inv(P)*transpose(B)*inv(B*inv(P)*transpose(B))*W;
X=dms2deg([69 59 50 ; 49 59 50 ; 60 00 20]);
A=[200 ; 500];
len_AB=500;
angle_AC=dms2deg([336 51 52]);
C=[A(1,1)+len_AB*sind(L(1,1))/sind(L(2,1))*sind(angle_AC) ; A(2,1)+len_AB*sind(L(1,1))/sind(L(2,1))*cosd(angle_AC)];
J=[500*(cosd(X(1,1))/sind(X(2,1))*sind(angle_AC)+sind(X(1,1))/sind(X(2,1))*cosd(angle_AC)) 500*(-sind(X(1,1))*cosd(X(2,1))/sind(X(2,1))/sind(X(2,1))*sind(angle_AC)+sind(X(1,1))/sind(X(2,1))*cosd(angle_AC)) ;
   500*(cosd(X(1,1))/sind(X(2,1))*cosd(angle_AC)-sind(X(1,1))/sind(X(2,1))*sind(angle_AC)) 500*(-sind(X(1,1))*cosd(X(2,1))/sind(X(2,1))/sind(X(2,1))*cosd(angle_AC)-sind(X(1,1))/sind(X(2,1))*sind(angle_AC))];
tmp=J*diag([(20/206265)^2 (20/206265)^2])*transpose(J);
SD=sqrt(tmp);

%Q5
B=[1 -1 0 ; 1 0 -1];
W=-[50.2/dms2rad([00 27 00])-120.5/dms2rad([01 05 00]) ; 50.2/dms2rad([00 27 00])-200.0/dms2rad([01 47 00])];
length=[50.2 ; 120.5 ; 200.0];
angle=[dms2rad([00 27 00]) ; dms2rad([01 05 00]) ; dms2rad([01 47 00])];
L=length./angle;
P=diag([1/((1/angle(1,1)*0.4)^2+(length(1,1)/angle(1,1)^2*dms2rad([00 0.2 00]))^2) 1/((1/angle(2,1)*0.4)^2+(length(2,1)/angle(2,1)^2*dms2rad([00 0.2 00]))^2) 1/((1/angle(3,1)*0.4)^2+(length(3,1)/angle(3,1)^2*dms2rad([00 0.2 00]))^2)]);
V=inv(P)*transpose(B)*inv(B*inv(P)*transpose(B))*W;
X=L+V;
Q=B*inv(P)*transpose(B);
tmp=transpose(V)*P*V/2*(inv(P)-inv(P)*transpose(B)*inv(Q)*B*inv(P));
SD=sqrt(tmp)

%Q22
B=[-1 -1 0 1 0 0 0 ; 0 -1 -1 0 1 0 0 ; -1 -1 -1 0 0 1 0 ; 1 1 1 0 0 0 1];
L_rad=[dms2rad([32 30 30]) ; dms2rad([20 01 20]) ; dms2rad([45 48 25]) ; dms2rad([52 32 30]) ; dms2rad([65 49 55]) ; dms2rad([98 20 40]) ; dms2rad([261 39 05])];
L_dms=rad2dms(L_rad);
W_rad=[L_rad(1,1)+L_rad(2,1)-L_rad(4,1) ; L_rad(2,1)+L_rad(3,1)-L_rad(5,1) ; L_rad(1,1)+L_rad(2,1)+L_rad(3,1)-L_rad(6,1) ; 2*pi()-L_rad(1,1)-L_rad(2,1)-L_rad(3,1)-L_rad(7,1)];
W_dms=rad2dms(W_rad);
P=diag([1 1 1 1 1 1 1]);
V_rad=inv(P)*transpose(B)*inv(B*inv(P)*transpose(B))*W_rad;
V_dms=rad2dms(V_rad);
X_rad=L_rad+V_rad;
X_dms=rad2dms(X_rad);
Q=B*inv(P)*transpose(B);
tmp=transpose(V_dms(:,end))*P*V_dms(:,end)/4*(inv(P)-inv(P)*transpose(B)*inv(Q)*B*inv(P));
SD=sqrt(tmp);

%Q23
L=[10.008 ; 19.987 ; 10.012 ; -15.004 ; -4.992 ; 10.010 ; 4.998];
B=[1 1 0 1 0 0 1 ; 0 1 -1 0 0 -1 0 ; 0 0 0 1 -1 1 0];
W=[20-L(1,1)-L(2,1)-L(4,1)-L(7,1) ; -L(2,1)+L(3,1)+L(6,1) ; -L(4,1)+L(5,1)-L(6,1)];
V=B.'*inv(B*B.')*W;
X=L+V;
H=[50+X(1,1)+X(2,1) ; 50+X(1,1)+X(2,1)+X(4,1) ; 50+X(1,1)+X(3,1) ; 50+X(1,1)];
Q=B*B.';
tmp=V.'*V/3*(inv(eye(7))-B.'*inv(Q)*B);
SD_LL=sqrt(tmp);
SD_P=sqrt([1 1]*[tmp(1,1) tmp(1,2) ; tmp(1,2) tmp(2,2)]*[1 ; 1]);
SD_Q=sqrt([-1]*[tmp(7,7)]*[-1]);
SD_R=sqrt([1 1]*[tmp(1,1) tmp(1,3) ; tmp(1,3) tmp(3,3)]*[1 ; 1]);
SD_S=sqrt([1]*[tmp(1,1)]*[1]);
SD_PR=sqrt(tmp(6,6));
SD_QS=sqrt([1 1]*[tmp(3,3) tmp(3,5) ; tmp(5,3) tmp(5,5)]*[1 ; 1]);

%Q24
B=[1 -1 0 0 ; 1 0 -1 0 ; 1 0 0 -1];
h=[72.840 ; -16.772 ; 14.802 ; 39.774];
H=[110.394 ; 200.052 ; 168.478 ; 143.506];
W=[H(2,1)-H(1,1)-h(1,1)+h(2,1) ; H(3,1)-H(1,1)-h(1,1)+h(3,1) ; H(4,1)-H(1,1)-h(1,1)+h(4,1)];
P=diag([1/25 ; 1/18 ; 1/15 ; 1/12]);
V=inv(P)*transpose(B)*inv(B*inv(P)*transpose(B))*W;
X=H+h+V;
tmp=V.'*P*V/3;
sigma0=sqrt(tmp);
Q=B*inv(P)*transpose(B);
tmp=transpose(V)*P*V/3*(inv(P)-inv(P)*transpose(B)*inv(Q)*B*inv(P));
SD=sqrt(tmp);





