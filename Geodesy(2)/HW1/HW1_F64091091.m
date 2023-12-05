%%% Initial Value %%%
a = 6378137; % semimajor axis (m)
f = 1/298.257222101; % flattening 
b = a - a*f; % semiminor axis (m)

re = 9.7803267715; % normal gravity at equator (m*s^-2)
e1 = sqrt(a^2-b^2) / a; % first excentricity
w = 7.292115*10^-5; % angular velocity (rad*s^-1)
GM = 3986005*10^8; % geocentric gravitational canstant
k = 0.001931851353;
B = -90: 1: 90;

r1 = 0.0052790414;
r2 = 0.0000232718;
r3 = 0.0000001262;
r4 = 0.0000000007;

%%% The Normal Gravity on the Level Ellipsoid within the Grid of One Degree - 1  %%%
R1 = zeros(181, 361);
for x = 0: 1: 180
    for y = 0: 1: 360
        R1(x+1, y+1) = re * (1 + r1*(sind(B(x+1))^2) + r2*(sind(B(x+1))^4) + r3*(sind(B(x+1))^6) + r4*(sind(B(x+1))^8));
    end
end
title_ = ['The Normal Gravity on the Level Ellipsoid', sprintf('\n'), 'within the Grid of One Degree by Using Formula 1'];
get_2Dplot(R1, title_, 'longtitude', 'latitude', 'normal_gravity_ellipsoid_1.png');

%%% The Normal Gravity on the Level Ellipsoid within the Grid of One Degree - 2 %%%
R2 = zeros(181, 361);
for x = 0: 1: 180
    for y = 0: 1: 360
        R2(x+1, y+1) = re * ((1 + k*(sind(B(x+1))^2)) / sqrt((1 - (e1^2)*(sind(B(x+1))^2))));
    end
end
title_ = ['The Normal Gravity on the Level Ellipsoid', sprintf('\n'), 'within the Grid of One Degree by Using Formula 2'];
get_2Dplot(R2, title_, 'longtitude', 'latitude', 'normal_gravity_ellipsoid_2.png');

%%% The Difference between Two Formula of The Normal Gravity on the Level Ellipsoid within the Grid of One Degree %%%
difference = zeros(181, 361);
for x = 0: 1: 180
    for y = 0: 1: 360
        difference(x+1, y+1) = R1(x+1, y+1) - R2(x+1, y+1);
    end
end
title_ = ['The Difference between Two Formula of The Normal Gravity',  sprintf('\n'), 'on the Level Ellipsoid within the Grid of One Degree'];
get_2Dplot(difference, title_, 'longtitude', 'latitude', 'difference_ellipsoid.png');

%%% The Normal Gravity Potential on the Level Ellipsoid within the Grid of One Degree %%%
u = b;
U = zeros(181, 361);
for x = 0:180
    for y = 0:360
        U(x+1, y+1) = (GM/e1)*atan(e1/u) + ((w^2)*(a^2)) / 3;
    end
end
title_ = ['The Normal Gravity Potential on the Level',  sprintf('\n'), 'Ellipsoid within the Grid of One Degree'];
get_2Dplot(U, title_, 'longtitude', 'latitude', 'normal_gravity_potential_ellipsoid.png');

%%% The Normal Gravity above the Level Ellipsoidal Height 10 Meter within the Grid of One Degree - 1 %%%
R1 = zeros(181, 361);
m = (w^2) * (a^2) * b / GM;
h = 10;
for x = 0: 1: 180
    for y = 0: 1: 360
        R1(x+1, y+1) = re * (1 + r1*(sind(B(x+1))^2) + r2*(sind(B(x+1))^4) + r3*(sind(B(x+1))^6) + r4*(sind(B(x+1))^8));
    end
end
rh1 = R1 * (1 - (2/a)*(1+f+m-2*f*(sind(x+1)^2)*h + 3*h^2/(a+2)));
title_ = ['The Normal Gravity above the Level Ellipsoidal', sprintf('\n'), 'Height 10 Meter within the Grid of One Degree by Using Formula 1'];
get_2Dplot(rh1, title_, 'longtitude', 'latitude', 'normal_gravity_height_1.png');

%%% The Normal Gravity above the Level Ellipsoidal Height 10 Meter within the Grid of One Degree - 2 %%%
R2 = zeros(181, 361);
m = (w^2) * (a^2) * b / GM;
h = 10;
for x = 0: 1: 180
    for y = 0: 1: 360
        R2(x+1, y+1) = re * ((1 + k*(sind(B(x+1))^2)) / sqrt((1 - (e1^2)*(sind(B(x+1))^2))));
    end
end
rh2 = R2 * (1 - (2/a)*(1+f+m-2*f*(sind(x+1)^2)*h + 3*h^2/(a+2)));
title_ = ['The Normal Gravity above the Level Ellipsoidal', sprintf('\n'), 'Height 10 Meter within the Grid of One Degree by Using Formula 2'];
get_2Dplot(rh2, title_, 'longtitude', 'latitude', 'normal_gravity_height_2.png');

%%% The Difference between Two Formula of The Normal Gravity above the Level Ellipsoidal Height 10 Meter within the Grid of One Degree %%%
difference = zeros(181, 361);
for x = 0: 1: 180
    for y = 0: 1: 360
        difference(x+1, y+1) = R1(x+1, y+1) - R2(x+1, y+1);
    end
end
title_ = ['The Difference between Two Formula of The Normal Gravity above the',  sprintf('\n'), 'Level Ellipsoidal Height 10 Meter within the Grid of One Degree'];
get_2Dplot(difference, title_, 'longtitude', 'latitude', 'difference_height.png');

%%% The Normal Gravity Potential above the Level Ellipsoidal Height 10 Meter within the Grid of One Degree %%%
u = b + 10;
q = (1/2) * ((1+3*((u^2)/(e1^2)))*atan(e1/u) - 3*(u/e1));
q0 = (1/2) * ((1+3*((b^2)/(e1^2)))*atan(e1/b) - 3*(b/e1));
for x = 0: 1: 180
    for y = 0: 1 :360
        U(x+1, y+1) = (GM/e1)*atan(e1/u) + ((w^2)/2)*(a^2)*(q/q0)*((sind(B(x+1)))^2-1/3) + ((w^2)/2)*(u^2+e1^2)*((cosd(B(x+1)))^2);
    end
end
title_ = ['The Normal Gravity Potential above the Level', sprintf('\n'), 'Ellipsoidal Height 10 Meter within the Grid of One Degree'];
get_2Dplot(U, title_, 'longtitude', 'latitude', 'normal_gravity_potential_height.png');

%%% Function %%%
function output = get_2Dplot(z_, title_, xlabel_, ylabel_, filename_)
    imagesc(z_);
    set(gca, 'XTick', (1: 45: 361), 'XTicklabel', (0: 45: 360));
    set(gca, 'YTick', (1: 10: 181), 'YTicklabel', (90: -10: -90));
    title(title_);
    xlabel(xlabel_);
    ylabel(ylabel_);
    colorbar;

    saveas(gcf, filename_);
    output = [];
end