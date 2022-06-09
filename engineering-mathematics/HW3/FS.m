function [an , bn , f] = fourier_series(fx , x , n , a , b)
    
    clc;
    clf;
    clear all;

    if nargin == 3
        a = -pi();
        b = pi();
    end

    l = (b-a)/2;

    if a+b
        fx = subs(fx , x , x+1+a);
    end
    
    an = int(fx , x , -1 , 1)/1;
    bn = [];
    f = an/2;

    for i = 1:n
        ann = int(fx*cos(i*pi()*x/1) , x , -1 , 1)/1;
        bnn = int(fx*sin(i*pi()*x/1) , x , -1 , 1)/1;
        an = [an , ann];
        bn = [bn , bnn];
        f = f+ann*coscos(i*pi()*x/1)+bnn*sin(i*pi()*x/1);
    end
    
    if a+b
        f = subs(fx , x , x-1-a);
    end    
    
end




