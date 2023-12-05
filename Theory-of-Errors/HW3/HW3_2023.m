clc;
clear;
%%% initial value %%%
syms HC HD HE;

h1 = 23.9900;
h2 = 11.9946;
h3 = 11.9933;
h4 = 06.0008;
h5 = 17.9993;
h6 = 06.0011;
h7 = 03.9988;
h8 = 10.0004;

HA = 100.0000;
HB = 108.0000;
file = fopen('HW3_output.txt', 'w');

%%% mixed model - 1 %%%
n = 8;
u = 4;
syms dHC dHD dHE da HC HD HE a;
HC0 = 118.0004;
HD0 = 123.9900;
HE0 = 111.9988;
a0 = 1.0000;

fun_A = [ 000 - dHD + 000 + h1*da;
          000 - dHD + dHE + h2*da;
          000 + 000 - dHE + h3*da;
          dHC - dHD + 000 + 00000;
         -dHC + 000 + 000 + 00000;
         -dHC + 000 + dHE + 00000;
          000 + 000 - dHE + 00000;
         -dHC + 000 + 000 + 00000];

A = zeros(size(fun_A, 1), u);
for i = 1: size(fun_A, 1)
    A(i, u-3) = diff(fun_A(i), dHC, 1);
    A(i, u-2) = diff(fun_A(i), dHD, 1);
    A(i, u-1) = diff(fun_A(i), dHE, 1);
    A(i, u) = diff(fun_A(i), da, 1);
end
B = diag([a, a, a, 1, 1, 1, 1, 1]);
W = [HD - HA - a*h1;
     HD - HE - a*h2;
     HE - HA - a*h3;
     HD - HC - h4;
     HC - HA - h5;
     HC - HE - h6;
     HE - HB - h7;
     HC - HB - h8];

flag = 0;
while true
     if flag == 0
          X = [HC0; HD0; HE0; a0];
          flag = 1;
     end
     W_ = subs(W, [HC; HD; HE; a], X);
     B_ = subs(B, a, X(4));
     dx = vpa(inv(transpose(A)*inv(B_*transpose(B_))*A)*transpose(A)*inv(B_*transpose(B_))*W_);
     X = X + dx;
     if max(abs(double(dx)) < 10^(-6))
          break;
     end
end
X_hat = X;

v_hat = transpose(B_)*inv(B_*transpose(B_))*(W_-A*dx);
sigma0_hat = sqrt((transpose(v_hat)*v_hat)/(n-u));
sigmaHH_hat = sigma0_hat^2 * inv(transpose(A)*A);
sigmaHC_hat_mixed1 = sqrt(sigmaHH_hat(1, 1));
sigmaHD_hat_mixed1 = sqrt(sigmaHH_hat(2, 2));
sigmaHE_hat_mixed1 = sqrt(sigmaHH_hat(3, 3));
sigmaa_hat = sqrt(sigmaHH_hat(4, 4));

fprintf(file, '%% ----- Mixed Model - 1 ----- %%\n');
fprintf(file, 'HC_hat = %7.4f (m), sigma_HC = %5.4f (m)\n', X_hat(1), sigmaHC_hat_mixed1);
fprintf(file, 'HD_hat = %7.4f (m), sigma_HD = %5.4f (m)\n', X_hat(2), sigmaHD_hat_mixed1);
fprintf(file, 'HE_hat = %7.4f (m), sigma_HE = %5.4f (m)\n', X_hat(3), sigmaHE_hat_mixed1);
fprintf(file, 'a_hat = %6.5f, sigma_a = %7.6f\n\n', X_hat(4), sigmaa_hat);

%% GMM model - 1 %%
n = 8;
u = 3;

L = [h1 + HA + 00;
     h2 + 00 + 00;
     h3 + HA + 00;
     h4 + 00 + 00;
     h5 + HA + 00;
     h6 + 00 + 00;
     h7 + 00 + HB;
     h8 + 00 + HB];

fun = [ 00 + HD + 00;
        00 + HD - HE;
        00 + 00 + HE;
       -HC + HD + 00; 
        HC + 00 + 00; 
        HC + 00 - HE; 
        00 + 00 + HE; 
        HC + 00 + 00];
A = zeros(size(fun, 1), u);
for i = 1: size(fun, 1)
    A(i, u-2) = diff(fun(i), HC, 1);
    A(i, u-1) = diff(fun(i), HD, 1);
    A(i, u) = diff(fun(i), HE, 1);
end
X_hat = inv(transpose(A)*A)*transpose(A)*L;

v_hat = A*X_hat - L;
sigma0_hat = sqrt((transpose(v_hat)*v_hat)/(n-u));
sigmaHH_hat = sigma0_hat^2 * inv(transpose(A)*A);
sigmaHC_hat_GMM1 = sqrt(sigmaHH_hat(1, 1));
sigmaHD_hat_GMM1 = sqrt(sigmaHH_hat(2, 2));
sigmaHE_hat_GMM1 = sqrt(sigmaHH_hat(3, 3));

precision_update = [1-sigmaHC_hat_mixed1/sigmaHC_hat_GMM1; 1-sigmaHD_hat_mixed1/sigmaHD_hat_GMM1; 1-sigmaHE_hat_mixed1/sigmaHE_hat_GMM1];

fprintf(file, '%% ----- GMM Model - 1 ----- %%\n');
fprintf(file, 'HC_hat = %7.4f (m), sigma_HC = %5.4f (m)\n', X_hat(1), sigmaHC_hat_GMM1);
fprintf(file, 'HD_hat = %7.4f (m), sigma_HD = %5.4f (m)\n', X_hat(2), sigmaHD_hat_GMM1);
fprintf(file, 'HE_hat = %7.4f (m), sigma_HE = %5.4f (m)\n\n', X_hat(3), sigmaHE_hat_GMM1);
fprintf(file, 'HC_hat_update = %3.1f%%\n', precision_update(1)*100);
fprintf(file, 'HD_hat_update = %3.1f%%\n', precision_update(2)*100);
fprintf(file, 'HE_hat_update = %3.1f%%\n\n', precision_update(3)*100);

%%% mixed model - 2 %%%
s = 8;
u = 3;
p = 1;
syms dHC dHD dHE da HC HD HE a;

fun_A = [ 00 + HD + 00;
          00 + HD - HE;
          00 + 00 + HE;
         -HC + HD + 00; 
          HC + 00 + 00; 
          HC + 00 - HE; 
          00 + 00 + HE; 
          HC + 00 + 00];

A = zeros(size(fun_A, 1), u);
for i = 1: size(fun_A, 1)
    A(i, u-2) = diff(fun(i), HC, 1);
    A(i, u-1) = diff(fun(i), HD, 1);
    A(i, u) = diff(fun(i), HE, 1);
end

B = diag([-1, -1, -1, -1, -1, -1, -1, -1]);
W = [h1 + HA;
     h2 + 00;
     h3 + HA;
     h4 + 00;
     h5 + HA;
     h6 + 00;
     h7 + HB;
     h8 + HB];

H = [0, 1 ,0];
h = 15.9995 + HB;

Qww = B*transpose(B);
Qyy = transpose(A)*inv(Qww)*A;
Qzz = H*inv(Qyy)*transpose(H);
X_hat = inv(Qyy)*transpose(A)*inv(Qww)*W + inv(Qyy)*transpose(H)*inv(H*inv(Qyy)*transpose(H))*(h-H*inv(Qyy)*transpose(A)*inv(Qww)*W)

ka_hat = inv(Qww)*(W-A*X_hat);
v_hat = transpose(B)*ka_hat
sigma0_hat = sqrt((transpose(v_hat)*v_hat)/(s+p-u))
sigmaHH_hat = sigma0_hat^2*inv(Qyy) - sigma0_hat^2*inv(Qyy)*transpose(H)*inv(Qzz)*H*inv(Qyy)
sigmaHC_hat_mixed2 = sqrt(sigmaHH_hat(1, 1))
sigmaHD_hat_mixed2 = sqrt(sigmaHH_hat(2, 2))
sigmaHE_hat_mixed2 = sqrt(sigmaHH_hat(3, 3))

fprintf(file, '%% ----- Mixed Model - 2 ----- %%\n');
fprintf(file, 'HC_hat = %7.4f (m), sigma_HC = %5.4f (m)\n', X_hat(1), sigmaHC_hat_mixed2);
fprintf(file, 'HD_hat = %7.4f (m), sigma_HD = %5.4f (m)\n', X_hat(2), sigmaHD_hat_mixed2);
fprintf(file, 'HE_hat = %7.4f (m), sigma_HE = %5.4f (m)\n\n', X_hat(3), sigmaHE_hat_mixed2);

%% GMM model - 2 %%
n = 8;
u = 2;

HD = 15.9995 + HB;
L = [h1 + HA + 00 - HD;
     h2 + 00 + 00 - HD;
     h3 + HA + 00 + 00;
     h4 + 00 + 00 - HD;
     h5 + HA + 00 + 00;
     h6 + 00 + 00 + 00;
     h7 + 00 + HB + 00;
     h8 + 00 + HB + 00];

fun = [ 00 + 00;
        00 - HE;
        00 + HE; 
       -HC + 00; 
        HC + 00; 
        HC - HE; 
        00 + HE;  
        HC + 00];

A = zeros(size(fun, 1), u);
for i = 1: size(fun, 1)
    A(i, u-1) = diff(fun(i), HC, 1);
    A(i, u) = diff(fun(i), HE, 1);
end
X_hat = inv(transpose(A)*A)*transpose(A)*L;

v_hat = A*X_hat - L;
sigma0_hat = sqrt((transpose(v_hat)*v_hat)/(n-u));
sigmaHH_hat = sigma0_hat^2 * inv(transpose(A)*A);
sigmaHC_hat_GMM2 = sqrt(sigmaHH_hat(1, 1));
sigmaHE_hat_GMM2 = sqrt(sigmaHH_hat(2, 2));

fprintf(file, '%% ----- GMM Model - 2 ----- %%\n');
fprintf(file, 'HC_hat = %7.4f (m), sigma_HC = %5.4f (m)\n', X_hat(1), sigmaHC_hat_GMM2);
fprintf(file, 'HD_hat = %7.4f (m), sigma_HD = %5.4f (m)\n', HD, 0);
fprintf(file, 'HE_hat = %7.4f (m), sigma_HE = %5.4f (m)\n', X_hat(2), sigmaHE_hat_GMM2);