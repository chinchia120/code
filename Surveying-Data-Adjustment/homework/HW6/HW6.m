%clc;
%clear all;
format long;

%Q19.12
syms XB YB dXB dYB;

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

fi_AMk1 = atan2(XMk1-XA,YMk1-YA);
fi_AB = fi_AMk1+ThetaA-2*pi();
XB0 = XA+LAB*sin(fi_AB);
YB0 = YA+LAB*cos(fi_AB);

F1 = LAB+((XB-XA)/LAB)*dXB+((YB-YA)/LAB)*dYB;
F2 = LBC-((XC-XB)/LBC)*dXB-((YC-YB)/LBC)*dYB;
F3 = ThetaA+((YB-YA)/LAB^2)*dXB-((XB-XA)/LAB^2)*dYB;
F4 = ThetaB+(-(YC-YB)/LBC^2+(YA-YB)/LAB^2)*dXB+((XC-XB)/LBC^2-(XA-XB)/LAB^2)*dYB;
F5 = ThetaC-((YB-YC)/LBC^2)*dXB+((XB-XC)/LBC^2)*dYB;

A = [diff(F1,dXB,1) diff(F1,dYB,1) ; 
    diff(F2,dXB,1) diff(F2,dYB,1) ; 
    diff(F3,dXB,1) diff(F3,dYB,1) ; 
    diff(F4,dXB,1) diff(F4,dYB,1) ; 
    diff(F5,dXB,1) diff(F5,dYB,1)];
L = [LAB ; LBC ; ThetaA ; ThetaB ; ThetaC];    
X = [sqrt((XB-XA)^2+(YB-YA)^2) ; 
    sqrt((XB-XC)^2+(YB-YC)^2) ; 
    2*pi()+(atan2(XB-XA,YB-YA)-atan2(XMk1-XA,YMk1-YA)) ; 
    atan2(XC-XB,YC-YB)-atan2(XA-XB,YA-YB) ; 
    atan2(XMk2-XC,YMk2-YC)-atan2(XB-XC,YB-YC)];
W = L-X;
P = diag([1/SLAB^2 1/SLBC^2 1/SThetaA^2 1/SThetaB^2 1/SThetaC^2]);

i = 1;
while 1
    if i == 1
        XB_ = XB0;
        YB_ = YB0;
    end
    dC = vpa(subs(inv(A.'*P*A)*A.'*P*W,[XB YB],[XB_ YB_]));
    XB_ = XB_+dC(1,1);
    YB_ = YB_+dC(2,1);
    i = i+1;
    if dC(1,1)<10^-6 && dC(2,1)<10^-6
        break;
    end    
end    

X_ = vpa(subs(X,[XB YB],[XB_ YB_]));
V = vpa(X_-L);
SD_0 = vpa(sqrt((V.'*P*V/(5-2)))); %ans_Q1

XB_adj = XB_; %ans_Q2
YB_adj = YB_; %ans_Q2
A_ = subs(A,[XB YB],[XB_ YB_]);
Cov_B = SD_0^2*inv(A_.'*P*A_);
SD_XB = vpa(sqrt(Cov_B(1,1))); %ans_Q2
SD_YB = vpa(sqrt(Cov_B(2,2))); %ans_Q2

X_adj = vpa(X_); %ans_Q3
V_adj = vpa(X_-L); %ans_Q3
J = [diff(X(1,1),XB,1) diff(X(1,1),YB,1) ; 
    diff(X(2,1),XB,1) diff(X(2,1),YB,1) ;
    diff(X(3,1),XB,1) diff(X(3,1),YB,1) ;
    diff(X(4,1),XB,1) diff(X(4,1),YB,1) ;
    diff(X(5,1),XB,1) diff(X(5,1),YB,1)];
Cov_obs = subs(J*Cov_B*J.',[XB YB],[XB_ YB_]);
SD_LAB = vpa(sqrt(Cov_obs(1,1))); %ans_Q3
SD_LBC = vpa(sqrt(Cov_obs(2,2))); %ans_Q3
SD_ThetaA = vpa(sqrt(Cov_obs(3,3))); %ans_Q3
SD_ThetaB = vpa(sqrt(Cov_obs(4,4))); %ans_Q3
SD_thetaC = vpa(sqrt(Cov_obs(5,5))); %ans_Q3

N_inv = vpa(inv(A_.'*P*A_)); %ans_Q4

X_0025_3 = 9.348;
X_0975_3 = 0.216;
X2_test = 3*SD_0^2/1;
while 1
    if X2_test>X_0975_3 && X2_test<X_0025_3 %ans_Q5
        fprintf('Can not reject\n');
        break;
    end
    fprintf('Reject\n');    
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

i = 0;
while 1
    if i == 0
        a_ = a0;
        b_ = b0;
    end
    dX = vpa(subs(inv(A.'*A)*A.'*W,[a b],[a_ b_]));
    a_ = a_ + dX(1,1);
    b_ = b_ + dX(2,1);
    i = 1;
    if dX(1,1)<10^-6 && dX(2,1)<10^-6
        break;
    end
end

F1_ = vpa(subs(F1,a,a_));
F2_ = vpa(subs(F2,b,b_));
F3_ = vpa(subs(F3,[a b],[a_ b_]));
F4_ = vpa(subs(F4,[a b],[a_ b_]));
F5_ = vpa(subs(F5,[a b],[a_ b_]));

dms_a = vpa(rad2dms(F1_));
dms_b = vpa(rad2dms(F2_));
dms_c = vpa(rad2dms(F3_));

fi_AB_ = atan2(EB-EA,NB-NA);
fi_AC = fi_AB_+(2*pi()-F3);
fi_AC_ = vpa(fi_AB_+(2*pi()-F3_));

dms_fi_AB = vpa(rad2dms(fi_AB_));
dms_fi_AC = vpa(rad2dms(fi_AC_));

C = [EA+F4_*sin(fi_AC_) NA+F4_*cos(fi_AC_)];

A_ = vpa(subs(A,[a,b],[a_ b_]));
Cov_ab = (20/206265)^2*inv(A_.'*A_);

EC = EA+F4*sin(fi_AC);
NC = NA+F4*cos(fi_AC);
J = [diff(EC,a,1) diff(EC,b,1) ; diff(NC,a,1) diff(NC,b,1)];
J_ = vpa(subs(J,[a b],[a_ b_]));
Cov_C = J_*Cov_ab*J_.';
SD_C2 = sqrt(Cov_C);

pre_EC = 100-SD_C2(1,1)/0.053*100;
pre_NC = 100-SD_C2(2,2)/0.047*100;

I_max = 0.5*(Cov_C(1,1)+Cov_C(2,2)+sqrt((Cov_C(1,1)-Cov_C(2,2))^2+4*Cov_C(1,2)^2));
I_min = 0.5*(Cov_C(1,1)+Cov_C(2,2)-sqrt((Cov_C(1,1)-Cov_C(2,2))^2+4*Cov_C(1,2)^2));

u = sqrt(I_max*5.991);
v = sqrt(I_min*5.991);
t = 0.5*atan2(-2*Cov_C(1,2),Cov_C(1,1)-Cov_C(2,2));
deg_t = vpa(rad2deg(t));

syms EC NC;

D = [100 -200];

%CD_ = pdist([C ; D],'euclidean');
CD = sqrt((D(1,1)-EC)^2+(D(1,2)-NC)^2);
J = [diff(CD,EC,1) diff(CD,NC,1)];
Cov_CD = vpa(subs(J*Cov_C*J.',[EC NC],[C(1,1) C(1,2)]));
SD_CD = sqrt(Cov_CD);

fi_CD = atan((D(1,1)-EC)/(D(1,2)-NC));
fi_CD_ = vpa(subs(fi_CD,[EC NC],[C(1,1) C(1,2)]));
dms_CD = vpa(rad2dms(fi_CD_));
J = [diff(fi_CD,EC,1) diff(fi_CD,NC,1)];
Cov_fi_CD = vpa(subs(J*Cov_C*J.',[EC NC],[C(1,1) C(1,2)]));
SD_fi_CD = sqrt(Cov_fi_CD);