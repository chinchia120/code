clc;
clear;

%Q4
pi_ = pi();
syms x n pi;
fx = x^2;
a0 = 1/(2*pi)*int(fx,x,-pi,pi);
an = 1/(pi)*int(fx*cos(n*x),x,-pi,pi)
bn = 1/(pi)*int(fx*sin(n*x),x,-pi,pi);

for i = 1:10
    subs(an,[n,pi],[i,pi_])
end    