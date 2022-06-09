clc;
clf;
clear;

syms x;
fx = x*(x-pi())*(x-2*pi());
left = 2;
size = 0.001;
right = 4;
width = right - left;

%T = 1
num_period = 1;
t = [left:size:right-size];
z = subs(fx,x,t);
array0 = [];

for i = 1:num_period
    array0 = [array0 z];
end

t = [left-(floor(num_period/2))*width:size:right+(floor(num_period/2))*width-size];
subplot(5,1,1);
plot(t,array0);

%T = 3
num_period = 3;
t = [left:size:right-size];
z = subs(fx,x,t);
array0 = [];

for i = 1:num_period
    array0 = [array0 z];
end

t = [left-(floor(num_period/2))*width:size:right+(floor(num_period/2))*width-size];
subplot(5,1,2);
plot(t,array0);

%T = 5
num_period = 5;
t = [left:size:right-size];
z = subs(fx,x,t);
array0 = [];

for i = 1:num_period
    array0 = [array0 z];
end

t = [left-(floor(num_period/2))*width:size:right+(floor(num_period/2))*width-size];
subplot(5,1,3);
plot(t,array0);

%T = 7
num_period = 7;
t = [left:size:right-size];
z = subs(fx,x,t);
array0 = [];

for i = 1:num_period
    array0 = [array0 z];
end

t = [left-(floor(num_period/2))*width:size:right+(floor(num_period/2))*width-size];
subplot(5,1,4);
plot(t,array0);

%T = 9
num_period = 9;
t = [left:size:right-size];
z = subs(fx,x,t);
array0 = [];

for i = 1:num_period
    array0 = [array0 z];
end

t = [left-(floor(num_period/2))*width:size:right+(floor(num_period/2))*width-size];
subplot(5,1,5);
plot(t,array0);