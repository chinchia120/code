format short

%%% HW2 Q7 %%%
syms A_ B_ C_ D_;

A_97= [2586564.018; 274040.686];
B_97 = [2561432.582; 271464.628];
C_97 = [2672140.277; 278091.433];

A_2010 = [2586357.074; 274869.216];
B_2010 = [2561226.156; 272292.830];
C_2010 = [2671934.535; 278920.676];

fun1 = A_*A_97(1, 1) - B_*A_97(2, 1) + C_;
fun2 = A_*B_97(1, 1) - B_*B_97(2, 1) + C_;
fun3 = A_*C_97(1, 1) - B_*C_97(2, 1) + C_;
fun4 = B_*A_97(1, 1) + A_*A_97(2, 1) + D_;
fun5 = B_*B_97(1, 1) + A_*B_97(2, 1) + D_;
fun6 = B_*C_97(1, 1) + A_*C_97(2, 1) + D_;

A = [diff(fun1, A_, 1) diff(fun1, B_, 1) diff(fun1, C_, 1) diff(fun1, D_, 1);
    diff(fun2, A_, 1) diff(fun2, B_, 1) diff(fun2, C_, 1) diff(fun2, D_, 1);
    diff(fun3, A_, 1) diff(fun3, B_, 1) diff(fun3, C_, 1) diff(fun3, D_, 1);
    diff(fun4, A_, 1) diff(fun4, B_, 1) diff(fun4, C_, 1) diff(fun4, D_, 1);
    diff(fun5, A_, 1) diff(fun5, B_, 1) diff(fun5, C_, 1) diff(fun5, D_, 1);
    diff(fun6, A_, 1) diff(fun6, B_, 1) diff(fun6, C_, 1) diff(fun6, D_, 1)];

X = [A_; B_; C_; D_];

W = [A_2010(1, 1) ; B_2010(1, 1); C_2010(1, 1); A_2010(2, 1); B_2010(2, 1); C_2010(2, 1)];

X_ = vpa(inv(A.'* A)* A.'* W)

%%% HW2 Q7 %%%
parameter_2010 = [1.6*10^(-3); 1.9*10^(-3); 2.4*10^(-3); -0.02*10^(-9); 0.0; 0.0; 0.0];
parameter_rate_2010 = [0.0*10^(-3); 0.0*10^(-3); -0.1*10^(-3); -0.03*10^(-9); 0.0; 0.0; 0.0];

parameter_2018 = parameter_2010 + (8)*parameter_rate_2010;

cor_2014_at2018 = [-2860754.858; 4961095.852; 2798783.722];
cor_rate_2014_at2018 = [0.010; 0.010; 0.010];

cor_2008_at2018 = vpa(cor_2014_at2018 + parameter_2018(1:3, 1) + diag([parameter_2018(4, 1), parameter_2018(4, 1), parameter_2018(4, 1)])*cor_2014_at2018)
cor_rate_2008_at2018 = vpa(cor_rate_2014_at2018 + parameter_2018(1:3, 1) + diag([parameter_2018(4, 1), parameter_2018(4, 1), parameter_2018(4, 1)])*cor_rate_2014_at2018)

cor_2008_at2010 = cor_2008_at2018 + cor_rate_2008_at2018*(-8)
