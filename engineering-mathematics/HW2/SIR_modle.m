% Time unit: 1 h
beta = (7076506-7071382)/(7076506*5000*24);
gamma = 0.05/24;
dt = 0.1;               % 6 min
D = 123;                 % Simulate for D days
N_t = floor(D*24/dt);   % Corresponding no of hours

t = linspace(0, N_t*dt, N_t+1);
S = zeros(N_t+1, 1);
I = zeros(N_t+1, 1);
R = zeros(N_t+1, 1);

% Initial condition
S(1) = 7076506;
I(1) = 1;
R(1) = 0;

% Step equations forward in time
for n = 1:N_t
   S(n+1) = S(n) - dt*beta*S(n)*I(n); 
   I(n+1) = I(n) + dt*beta*S(n)*I(n) - dt*gamma*I(n);
   R(n+1) = R(n) + dt*gamma*I(n);
end

plot(t, S, t, I, t, R);
legend('S', 'I', 'R', 'Location','northwest');
xlabel('hours');