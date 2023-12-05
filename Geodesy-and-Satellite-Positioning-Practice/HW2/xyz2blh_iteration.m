function y = xyz2blh(X, Y, Z, a, f)

X = -3039190.822;
Y = 5009755.660;
Z = 2516244.457;
a = 6378137.0;
f = 1/298.257223563;

dd = 180/pi;
landa = atan2(Y,X);

e = sqrt(2*f-f^2);

%initial values of fi and h.
fi0 = acos(sqrt(X^2+Y^2)/sqrt(X^2+Y^2+Z^2));
N0  = a/sqrt(1-(e*sin(fi0))^2);
h0 = 0;

dfi = 1;
while (dfi > 1e-12)
   N0  = a/sqrt(1-(e*sin(fi0))^2);
   fi  = atan((Z + e^2*N0*sin(fi0))/sqrt(X^2+Y^2));
   dfi = abs(fi - fi0);
   fi0 = fi;
end
N   = a / sqrt(1-e^2*(sin(fi))^2);

if fi > pi/4 | fi < -pi/4
   h = sqrt(X^2 + Y^2)/cos(fi) -N;
else
   h = Z / sin(fi) - N + e^2 * N;
end
y = [fi*dd; landa*dd; h];
