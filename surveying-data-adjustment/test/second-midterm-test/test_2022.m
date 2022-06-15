clc;
clear;
format long;

syms  X Y dX dY;

SAB = 40.050;
SBC = 30.050;
SAC = 50.090;

X0 = 40;
Y0 = 30;

F1 = X+dX;
F2 = Y+dY;
F3 = sqrt(X^2+Y^2)+X/sqrt(X^2+Y^2)*dX+Y/sqrt(X^2+Y^2)*dY;

A = [diff(F1,dX,1) diff(F1,dY,1) ; diff(F2,dX,1) diff(F2,dY,1) ; diff(F3,dX,1) diff(F3,dY,1)];
L = [SAB ; SBC ; SAC];
x = [X ; Y ; sqrt(X^2+Y^2)];
W = L-x;
P = diag([1/0.01^2 1/0.01^2 1/0.01^2]);

i = 1
while 1
    if i == 1
        X_ = X0;
        Y_ = Y0;
    end

    d = vpa(subs(inv(A.'*P*A)*A.'*P*W,[X Y],[X_ Y_]));
    X_ = X_+d(1,1);
    Y_ = Y_+d(2,1);   

    i = i+1;
    if i == 3
        break;
    end
end

x_ = vpa(subs(x,[X Y],[X_ Y_]));
V = vpa(x_-L);
SD_0 = vpa(sqrt((V.'*P*V/(3-2))));
A_ = subs(A,[X Y],[X_ Y_]);
Cov_x = SD_0^2*inv(A_.'*P*A_);
SD_X = vpa(sqrt(Cov_x(1,1)));
SD_Y = vpa(sqrt(Cov_x(2,2)));