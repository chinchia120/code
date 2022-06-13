clc;
clf;
clear;

t = 0:1/100:4;
y = cos(2*pi*t)+cos(5*pi*t);
noise = randn(1,size(y,2));
y_noise = y+noise;
FY_noise = fft(y_noise);
FY = FY_noise;

for i = 1:size(FY_noise,2)
    if abs(FY(i))<100
        FY(i) = 0;
    end
end

iFY = ifft(FY);

subplot(4,1,1);
plot(t,y_noise);

subplot(4,1,2);
plot(t,FY_noise);

subplot(4,1,3);
plot(t,FY);

subplot(4,1,4);
plot(t,iFY);