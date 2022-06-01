clc;
clear all;
format long;

%Q19.12
%syms EA NA EB NB EC NC EMk1 NMk1 EMk2 NMk2 LAB LBC ThetaA ThetaB ThetaC SLAB SLBC SThetaA SThetaB SThetaC;

syms dXB dYB;

XA = 5572.32;
YA = 6208.30;
XC = 9552.58;
YC = 6349.45;
XMk1 = 6212.18;
YMk1 = 4956.83;
XMk2 = 11547.42;
YMk2 = 6518.82;
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

fi_AMk1 = atan((XMk1-XA)/(YMk1-YA))+pi();
fi_AB = fi_AMk1+ThetaA-2*pi();
XB = XA+LAB*sin(fi_AB);
YB = YA+LAB*cos(fi_AB);


F1 = LAB+(XB-XA)/LAB*dXB+(YB-YA)/LAB*dYB;
F2 = LBC+(XC-XB)/LBC*dXB+(YC-YB)/LBC*dYB;
F3 = ThetaA+((YB-YA)/LAB^2)*dXB-((XB-XA)/LAB^2)*dYB;
F4 = ThetaB+(-(YC-YB)/LBC^2+(YA-YB)/LAB^2)*dXB+((XC-XB)/LBC^2-(XA-XB)/LAB^2)*dYB;
F5 = ThetaC-((YB-YC)/LBC^2)*dXB+((XB-XC)/LBC^2)*dYB;

A = [diff(F1,dXB,1) diff(F1,dYB,1) ; 
    diff(F2,dXB,1) diff(F2,dYB,1) ; 
    diff(F3,dXB,1) diff(F3,dYB,1) ; 
    diff(F4,dXB,1) diff(F4,dYB,1) ; 
    diff(F5,dXB,1) diff(F5,dYB,1)];
W = [LAB ; F2-LBC ; F3-ThetaA ; F4-ThetaB ; F5-ThetaC]
P = diag([1/SLAB^2 1/SLBC^2 1/SThetaA^2 1/SThetaB^2 1/SThetaC^2]); 

X = vpa(inv(A.'*P*A)*A.'*P*W);

for i = 1:2
    if i == 1
        dXB_ = 0;
        dYB_ = 0;
    end
    dC = vpa(subs(inv(A.'*P*A)*A.'*P*W,[dXB dYB],[dXB_ dYB_]));
    dXB_ = dXB_+dC(1,1);
    dYB_ = dYB_+dC(2,1);
end    



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
for i = 1:10
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

dms_a = vpa(rad2dms(F1_));
dms_b = vpa(rad2dms(F2_));
dms_c = vpa(rad2dms(F3_));

fi_AB = atan((EB-EA)/(NB-NA));
fi_AC = fi_AB+(2*pi()-F3);
fi_AC_ = vpa(fi_AB+(2*pi()-F3_));

dms_fi_AB = vpa(rad2dms(fi_AB));
dms_fi_AC = vpa(rad2dms(fi_AC_));

C = [EA+F4_*sin(fi_AC_) NA+F4_*cos(fi_AC_)];

A_ = vpa(subs(A,[a,b],[a_ b_]));
%V = vpa([F1_-a0 ; F2_-b0 ; F3_-c0 ; F4_-AC ; F5_-BC]);
Cov_ab = (20/206265)^2*inv(A_.'*A_);

EC = EA+F4*sin(fi_AC);
NC = NA+F4*cos(fi_AC);
J = [diff(EC,a,1) diff(NC,a,1) ; diff(EC,b,1) diff(NC,b,1)];
J_ = vpa(subs(J,[a b],[a_ b_]));
Cov_C = J_*Cov_ab*J_.';
SD_C2 = sqrt(Cov_C);

SD_C1 = sqrt([0.0021 -0.0004 ; -0.0004 0.0017]);
ans = SD_C2./SD_C1;
