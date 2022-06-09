clc;
clf;
clear;

syms x;
fx = x*(x-2)*(x-4);
left = 2;
size = 0.001;
right = 4;
width = right - left;

for i = 1:1:5
    num_period = 2*i-1;
    t = [left:size:right-size];
    z = subs(fx,x,t);
    array0 = [];

    for j = 1:num_period
        array0 = [array0 z];
    end

    t = [left-(floor(num_period/2))*width:size:right+(floor(num_period/2))*width-size];
    subplot(5,1,i);
    plot(t,array0);
end
