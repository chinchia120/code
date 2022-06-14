clc;
clf;
clear;

pi_ = pi();
syms x n pi;

x1 = linspace(0,pi_,10);
y1 = x1*0-3;
x2 = linspace(pi_,2*pi_,10);
y2 = x2*0+1;
x3 = linspace(2*pi_,3*pi_,10);
y3 = x3*0-3;
x4 = linspace(3*pi_,4*pi_,10);
y4 = x4*0+1;

x1_ = linspace(-pi_,0,10);
y1_ = x1_*0+1;
x2_ = linspace(-2*pi_,-pi_,10);
y2_ = x2_*0-3;
x3_ = linspace(-3*pi_,-2*pi_,10);
y3_ = x3_*0+1;
x4_ = linspace(-4*pi_,-3*pi_,10);
y4_ = x4_*0-3;

subplot(5,1,1);
plot([x4_ x3_ x2_ x1_ x1 x2 x3 x4],[y4_ y3_ y2_ y1_ y1 y2 y3 y4]);
legend('Original Function');

a0 = 1/(2*pi)*(int(-3,x,0,pi)+int(1,x,-pi,0));
an = 1/pi*(int(1*cos(n*x),x,-pi,0))+1/pi*(int(-3*cos(n*x),x,0,pi));
bn = 1/pi*(int(1*sin(n*x),x,-pi,0))+1/pi*(int(-3*sin(n*x),x,0,pi));
x_ = linspace(-4*pi_,4*pi_);

An = 0;
Bn = 0;
for i = 1:201
    
    an_ = subs(an*cos(n*x),[n pi],[i pi_]);
    An = An+an_;
    bn_ = subs(bn*sin(n*x),[n pi],[i pi_]);
    Bn = Bn+bn_;
    if i == 3
        subplot(5,1,2);
        plot(x_,subs(a0+An+Bn,x,x_));
        legend('n = 3');
    end

    if i == 7
        subplot(5,1,3);
        plot(x_,subs(a0+An+Bn,x,x_));
        legend('n = 7');
    end

    if i == 101
        subplot(5,1,4);
        plot(x_,subs(a0+An+Bn,x,x_));
        legend('n = 101');
    end
    
    if i == 201
        subplot(5,1,5);
        plot(x_,subs(a0+An+Bn,x,x_));
        legend('n = 201');
    end
    xlim([-4*pi_ 4*pi_]);
    ylim([-5 5]);
end 
