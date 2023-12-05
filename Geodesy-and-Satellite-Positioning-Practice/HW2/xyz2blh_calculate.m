X = -3039190.822;
Y = 5009755.660;
Z = 2516244.457;
a = 6378137.0;
f = 1/298.257223563; %% f = (a-b)/a

dd = 180/pi;
b = a - a*f;
p  = sqrt(X^2 + Y^2);
theta = atan((Z*a) / (p*b));
e1 = sqrt(a^2-b^2) / a;
e2 = sqrt(a^2-b^2) / b;

psi = atan((Z + (e2^2)*b*(sin(theta)^3)) / (p - (e1^2)*a*(cos(theta)^3)));
lambda = atan2(Y, X);
N = a / sqrt(1-e1^2*(sin(psi))^2);
h = (p / cos(psi)) - N;

y = [psi*dd; lambda*dd; h]