%%% Initial Value %%%
J = 1.08263*10^-3/sqrt(5); 
C = zeros(2,4);                 % spherical harmonic coefficients
S = zeros(2,4);                 % spherical harmonic coefficients
C(2,1) = -484165.14*10^-9 + J;  % C20
C(2,2) = -000000.21*10^-9;      % C21
C(2,3) = +002439.38*10^-9;      % C22
C(3,1) = +000957.16*10^-9;      % C30
C(3,2) = +002029.99*10^-9;      % C31
C(3,3) = +000904.79*10^-9;      % C32
C(3,4) = +000721.32*10^-9;      % C33
S(2,2) = +000001.38*10^-9;      % S21
S(2,3) = -001400.27*10^-9;      % S22
S(3,2) = +000248.20*10^-9;      % S31
S(3,3) = -000619.01*10^-9;      % S32
S(3,4) = +001414.35*10^-9;      % S33

R = 6378137;                    % semimajor axis (m)
GM = 3986005*10^8;              % geocentric gravitational canstant (m^3*s^-2)

%%% Gravity Anomaly %%%
gravity_anomaly = zeros(181, 361);
P = zeros(3,4);                 %Legendre function
for lat = 1: 181   
    P(2, 1:3) = legendre(2, cosd(lat-1)); 
    P(3, 1:4) = legendre(3, cosd(lat-1));
	for n = 2: 3   
        for m = 1: n+1 
            if m > 1            % fully normalized Legendre functions and coefficients
                P(n, m) = (-1)^(m+1) * P(n,m) * sqrt((2*(2*n+1) * factorial(n-m+1)) / factorial(n+m-1));
            else
                P(n, m) = P(n,m) * sqrt(((2*n+1) * factorial(n-m+1)) / factorial(n+m-1));
            end
                                % spherical harmonic
            gravity_anomaly(lat, :) = gravity_anomaly(lat, :) + P(n,m) * (C(n,m)*cosd((m-1)*(0:360)) + S(n,m)*sind((m-1)*(0:360))) * (n-1);
        end
	end
end
gravity_anomaly = gravity_anomaly * GM/R^2;

title_ = ['The Gravity Anomaly on the Level',  sprintf('\n'), 'Ellipsoid within the Grid of One Degree'];
get_2Dplot(gravity_anomaly, title_, 'longtitude', 'latitude', 'gravity_anomaly.png');

%%% Geoid Undulation %%%
geoid_undulation = zeros(181,361);
P = zeros(3,4);                 %Legendre function
for lat = 1:181   
    P(2, 1:3) = legendre(2, cosd(lat-1)); 
    P(3, 1:4) = legendre(3, cosd(lat-1));
	for n = 2: 3   
        for m = 1: n+1 
            if m>1              % fully normalized Legendre functions and coefficients
                P(n,m) = (-1)^(m+1) * P(n,m) * sqrt((2*(2*n+1) * factorial(n-m+1)) / factorial(n+m-1));
            else
                P(n,m) = P(n,m) * sqrt(((2*n+1) * factorial(n-m+1)) / factorial(n+m-1));
            end
                                % spherical harmonic
            geoid_undulation(lat, :) = geoid_undulation(lat, :) + P(n,m) * (C(n,m)*cosd((m-1)*(0:360)) + S(n,m)*sind((m-1)*(0:360)));
        end
	end
end
geoid_undulation = geoid_undulation * R;

title_ = ['The Geoid Undulation on the Level',  sprintf('\n'), 'Ellipsoid within the Grid of One Degree'];
get_2Dplot(geoid_undulation, title_, 'longtitude', 'latitude', 'geoid_undulation.png');

%%%  Free-air Anomaly  %%%
re = 9.7803267715;              % normal gravity at equator (m*s^-2)
r1 = 0.0052790414;
r2 = 0.0000232718;
r3 = 0.0000001262;
r4 = 0.0000000007;

fi = 20;                        % latitude (deg)
g = 9.84;                       % gravity (m*s^-2)
H = 10.000;                     % height (m)
r = re * (1 + r1*(sind(fi)^2) + r2*(sind(fi)^4) + r3*(sind(fi)^6) + r4*(sind(fi)^8));
delta_g0 = g + 0.3086*H*10^(-5) - r;

%%% Function %%%
function output = get_2Dplot(z_, title_, xlabel_, ylabel_, filename_)
    imagesc(z_);
    set(gca, 'XTick', (1: 45: 361), 'XTicklabel', (0: 45: 360));
    set(gca, 'YTick', (1: 30: 181), 'YTicklabel', (90: -30: -90));
    title(title_);
    xlabel(xlabel_);
    ylabel(ylabel_);
    colorbar;

    saveas(gcf, filename_);
    output = [];
end