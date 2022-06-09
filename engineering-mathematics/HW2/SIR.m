function SIR
clear all;
clc;    
clf;    

iterations = 1;  % Sets initial interation count to 1;
pausetime = 0.01;  % Shows solutions at each time step. 
runtime = 30;    % Duration time of simulation.
g = (7565.538326-7155.469044)/402;
b = (7071382-7064833)*7079349.538/(7071382*402);

% =============== Initial conditions for s、i and r ==========================
initials = 7071382/7079349.538;
initiali = 402/7079349.538;
initialr = (7565.538326-7155.469044)/7079349.538;

deq1 = @(t,x) [-b*x(1)*x(2); b*x(1)*x(2)-g*x(2); g*x(2)];
[t,sol] = ode45(deq1,[0 runtime],[initials initiali initialr]);

arraysize = size(t);  % Sets time array size for the for loop.

%============ Solutions are plotted at each time step =====================

for i = 1 : max(arraysize)  
    subplot(4,1,1)
    title('S')
    plot(t(i,1),sol(i,1),'b.','markersize',10,'MarkerFaceColor','b')
    legend('S')
    hold on;
    axis([0 30 0 1])

    subplot(4,1,2)
    title('I')
    plot(t(i,1),sol(i,2),'r.','markersize',10,'MarkerFaceColor','b')
    legend('I')
    hold on;
    axis([0 30 0 1])
    
    subplot(4,1,3)
    title('R')
    plot(t(i,1),sol(i,3),'g.','markersize',10,'MarkerFaceColor','b')
    legend('R')
    hold on;
    axis([0 30 0 1])
  
    subplot(4,1,4)
    title('SIR')
    hold on;
    plot(t(i,1),sol(i,1),'b.','markersize',10,'MarkerFaceColor','b')
    plot(t(i,1),sol(i,2),'r.','markersize',10,'MarkerFaceColor','b')
    plot(t(i,1),sol(i,3),'g.','markersize',10,'MarkerFaceColor','b')
    legend('S', 'I', 'R')
    axis([0 30 0 1])
    hold off;
    
    iterations = iterations + 1;   % Adds 1 to the iteration count. 
    pause(pausetime)
end