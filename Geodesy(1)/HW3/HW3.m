format short;

%%% Initial Value %%%
G = 6.673 * 10^(-11); % m^(3) * kg^(-1) * s^(-2)
omega = 7.292115 * 10^(-5); %rad * s^(-1)

%%% New Folder %%%
[status, msg, msgID] = mkdir('picture');

%%% Q1 %%%
R = 6370;
delta = 1: 60;

Area = double(zeros(1, 60));
for i = 1: 60
    Area(i) = R^2 * cosd(10) * delta(i)/60 * delta(i)/60;
end

get_2Dplot(delta, Area, 'Surface Area in Different Latitude and Longtitude', 'delta(min)', 'surface area(km^2)', './picture//Q1.png');

%%% Q3 %%%
cor1 = [0, 2000, 1000];
cor2 = [0, -1000, -1000];
cor3 = [0, -2000, -1000];
corP = [0, 0, 0];

m1 = 3000;
m2 = 2000;
m3 = 3000;

corC = (m1*cor1 + m2*cor2 + m3*cor3)/(m1 + m2 + m3);

gravitational_potential_P_center = get_gravitational_potential_cor(G, (m1+m2+m3), corC, corP);

gravitational_potential_P_m1 = get_gravitational_potential_cor(G, m1, cor1, corP);
gravitational_acceleration_P_m1 = get_gravitational_acceleration_cor(G, m1, cor1, corP);

gravitational_potential_P_m2 = get_gravitational_potential_cor(G, m2, cor2, corP);
gravitational_acceleration_P_m2 = get_gravitational_acceleration_cor(G, m2, cor2, corP);

gravitational_potential_P_m3 = get_gravitational_potential_cor(G, m3, cor3, corP);
gravitational_acceleration_P_m3 = get_gravitational_acceleration_cor(G, m3, cor3, corP);

%%% Q5 %%%
M = 10^(7); % kg
r = 1000*10^3; %% m
gravitational_potential = get_gravitational_potential_length(G, M, r);
gravitational_acceleration = get_gravitational_acceleration_length(G, M, r);

p = r * sind(30);
centrifugal_potential = get_centrifugal_potential(omega, p);
centrifugal_acceleration = get_centrifugal_acceleration(omega, p);

gravity_potrntial = gravitational_potential + centrifugal_potential;
gravity_acceleration = gravitational_acceleration + centrifugal_acceleration;

%%% Q6 %%%
syms s d
s_ = 20 * 10^3; % m
d_ = 10 * 10^3; % m
m = 1000; % kg
V = (G*m)/sqrt(s^2 + d^2);

V_ = double(vpa(get_gravitational_potential_length(G, m, subs(sqrt(s^2 + d^2), [s, d], [s_, d_]))));

V_diff_d1 = double(vpa(subs(diff(V, d, 1), [s, d], [s_, d_])));
V_diff_d2 = double(vpa(subs(diff(V, d, 2), [s, d], [s_, d_])));

%%% Q11 %%%
theta = 0: 0.01: pi()/2;
lambda = 0: 0.04: 2*pi();

[x, y] = meshgrid(theta, lambda);
P = legendre(4, cos(theta));
[P_, lambda_] = meshgrid(P(4, :), cos(3*lambda));
z = P_.*lambda_;

get_3Dplot(x, y, z, 'P_4_,_3(cos\theta)cos(3\lambda)', '\theta', '\lambda', 'z', './picture//Q11.png');

%%% Function %%%
function output = get_length(cor1, cor2)
    output = sqrt((cor1(1, 1)-cor2(1, 1))^2 + (cor1(1, 2)-cor2(1, 2))^2 + (cor1(1, 3)-cor2(1, 3))^2);
end

function output = get_gravitational_potential_cor(G, m1, cor1, cor2)
    l = get_length(cor1, cor2);
    output = (G*m1/l);
end

function output = get_gravitational_potential_length(G, M, r)
    output = (G*M/r);
end

function output = get_gravitational_acceleration_cor(G, m1, cor1, cor2)
    l = get_length(cor1, cor2);
    output = (G*m1/l^2);
end

function output = get_gravitational_acceleration_length(G, M, r)
    output = (G*M/r^2);
end

function output = get_centrifugal_potential(omega, p)
    output = (omega^2*p^2)/2;
end

function output = get_centrifugal_acceleration(omega, p)
    output = (omega^2*p);
end

function output = get_2Dplot(x, y, title_, xlabel_, ylabel_, filename_)
    plot(x, y);
    title(title_);
    xlabel(xlabel_);
    ylabel(ylabel_);

    saveas(gcf, filename_);
    output = [];
end

function output = get_3Dplot(x, y, z, title_, xlabel_, ylabel_, zlabel_, filename_)
    mesh(x, y, z)
    title(title_)
    xlabel(xlabel_)
    ylabel(ylabel_)
    zlabel(zlabel_)

    saveas(gcf, filename_)
    output = [];
end