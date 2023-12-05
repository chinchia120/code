%clc;
%clear all;
format long;

% Q6.17
syms pi;
syms W;
syms L;
syms s_W;
syms s_L;
J=[pi 2 ; (pi*W)/2+L W]
XX=[s_W^2 0 ; 0 s_L^2]
PA=J*XX*transpose(J)

% Q20
JJ=[-89.36690236 -519.6266275 ; 333.5218201 32.25071716]
X=[(20/206265)^2 0 ; 0 (20/206265)^2]
V=JJ*X*transpose(JJ)
 