clc;
clear all;
format long;

%Q19.12
%syms EA NA EB NB EC NC EMk1 NMk1 EMk2 NMk2 LAB LBC ThetaA ThetaB ThetaC SLAB SLBC SThetaA SThetaB SThetaC;

syms EB NB;

EA = 5572.32;
NA = 6208.30;
EC = 9552.58;
NC = 6349.45;
EMk1 = 6212.18;
NMk1 = 4956.83;
EMk2 = 11547.42;
NMk2 = 6518.82;
LAB = 2717.95;
LBC = 2589.28;
ThetaA = dms2rad([254 53 08]);
ThetaB = dms2rad([262 46 20]);
ThetaC = dms2rad([134 34 14]);
SLAB = 0.024;
SLBC = 0.024;
SThetaA = 5.3/206265;
SThetaB = 5.1/206265;
SThetaC = 5.2/206265;

F_LAB = sqrt((EB-EA)^2+(NB-NA^2));
F_LBC = sqrt((EC-EB)^2+(NC-NB^2));
F_ThetaA = atan((EB-EA)/(NB-NA))-atan((EMk1-EA)/(NMk1-NA));
F_ThetaB = atan((EC-EB)/(NC-NB))-atan((EA-EB)/(NA-NB));
F_ThetaC = atan((EMk2-EC)/(NMk2-NC))-atan((EB-EC)/(NB-NC));

A = [diff(F_LAB,EB,1) diff(F_LAB,NB,1) ; 
    diff(F_LBC,EB,1) diff(F_LBC,NB,1) ; 
    diff(F_ThetaA,EB,1) diff(F_ThetaA,NB,1) ; 
    diff(F_ThetaB,EB,1) diff(F_ThetaB,NB,1) ; 
    diff(F_ThetaC,EB,1) diff(F_ThetaC,NB,1)];
W = [F_LAB-LAB ; F_LBC-LBC ; F_ThetaA-ThetaA ; F_ThetaB-ThetaB ; F_ThetaC-ThetaC];
P = diag([1/SLAB^2 1/SLBC^2 1/SThetaA^2 1/SThetaB^2 1/SThetaC^2]); 

X = vpa(inv(A.'*P*A)*A.'*P*W)

%Q4
syms a b;

a0 = dms2rad([70 00 00]);
b0 = dms2rad([50 00 00]);
c0 = dms2rad([60 00 30]);
AB = 500;
AC = 613.353;
BC = 565.244;

EA = 200;
NA = 500;
EB = 500;
NB = 900;

F1 = a;
F2 = b;
F3 = pi()-a-b;
F4 = AB*sin(a)/sin(b);
F5 = AB*sin(a+b)/sin(b);

A = [1 0 ; 0 1 ; -1 -1 ; diff(F4,a,1) diff(F4,b,1) ; diff(F5,a,1) diff(F5,b,1)];
W = [0 ; 0 ; c0-a-b ; AC-F4 ; BC-F5];

for i = 1:3
    if i == 1
        a_ = a0;
        b_ = b0;
    end
    dX = vpa(subs(inv(A.'*A)*A.'*W,[a b],[a_ b_]));
    a_ = a_ + dX(1,1);
    b_ = b_ + dX(2,1);
end

F1_ = vpa(subs(F1,a,a_));
F2_ = vpa(subs(F2,b,b_));
F3_ = vpa(subs(F3,[a b],[a_ b_]));
F4_ = vpa(subs(F4,[a b],[a_ b_]));
F5_ = vpa(subs(F5,[a b],[a_ b_]));

fi_AB = atan((EB-EA)/(NB-NA));
fi_AC = fi_AB+(2*pi()-F3);
fi_AC_ = vpa(fi_AB+(2*pi()-F3_));
C = [EA+F4_*sin(fi_AC_) NA+F4_*cos(fi_AC_)];

A_ = vpa(subs(A,[a,b],[a_ b_]));
V = vpa([F1_-a0 ; F2_-b0 ; F3_-c0 ; F4_-AC ; F5_-BC]);
Cov_ab = ((V.'*V)/(5-2))^2*inv(A_.'*A_);

EC = EA+F4*sin(fi_AC);
NC = NA+F4*cos(fi_AC);
J = [diff(EC,a,1) diff(NC,a,1) ; diff(EC,b,1) diff(NC,b,1)];
J_ = vpa(subs(J,[a b],[a_ b_]));
Cov_C = J_*Cov_ab*J_.';
SD_C2 = sqrt(Cov_C);

SD_C1 = sqrt([0.0021 -0.0004 ; -0.0004 0.0017]);
ans = SD_C2./SD_C1;

git add -A
git commit -m "update"
git push