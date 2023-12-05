format short;

%%% New Folder %%%
[status, msg, msgID] = mkdir('picture');

%%% Q2 %%%
n = 91;
m = 90;
theta = 0: 0.01: pi()/2;
lambda = 0: 0.04: 2*pi();
[x, y] = meshgrid(theta, lambda);
P = legendre(n, cos(theta));
[P_, lambda_] = meshgrid(P(n, :), cos(m*lambda));
z = P_.*lambda_;

get_3Dplot(x, y, z, 'P_9_1_,_9_0(cos\theta)cos(90\lambda)', '\theta', '\lambda', 'z', './picture//Q2.png');

%%% Function %%%
function output = get_3Dplot(x, y, z, title_, xlabel_, ylabel_, zlabel_, filename_)
    mesh(x, y, z)
    title(title_)
    xlabel(xlabel_)
    ylabel(ylabel_)
    zlabel(zlabel_)

    saveas(gcf, filename_)
    output = [];
end