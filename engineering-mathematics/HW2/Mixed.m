% =========================================================================
% Matlab program to draw a phase portrait for the Mixed problem model. 
% The results are shown at each time step. 
% In addition, the user has the option of plotting a time series graph for x or y. 
% Equations are solved using a numerical ODE solver. 
%
%   
% Hone-Jay Chu 3/6/2020
% =========================================================================

function Mixed
clear  % Clears command window
clc    % Clears command history
clf    % Removes anything in the figure window before simulation. 

% ======================================================= 

countdown = 3;  % Countdown time length before simulation runs. 
iterations = 1;  % Sets initial interation count to 1;
pausetime = 0.5;  % Shows solutions at each time step. 
runtime = 120;    % Duration time of simulation.


% =============== Initial conditions for x and y ==========================
initialx1 = 0;
initialx2 = 150;

fprintf('----------------------------------\n4.1 Mixing problem model \n\nMatlab code written by James Adamss, modified by HJ Chu\n----------------------------------')
% fprintf('\n\nParameter values set,')
% %fprintf('\n\nalpha = %2.6f \nbeta = %2.6f \ngamma = %2.6f \ndelta = %2.6f ',alpha,beta,gamma,delta)
% fprintf('\n\nTo plot a time series graph for x, set the choice parameter in the code to 1.') 
% fprintf('\nSet choice to 2 to plot a time series graph for y.\n\nSimulation will run in ')

% for i = 3:-1:1
%     fprintf('\n%8i',countdown')    % Countdown to simulation start.
%     countdown = countdown-1;
%     pause(1)
% end

% Solves equations using numerical ODE solver 45 (nonstiff runge kutta)
deq1=@(t,x) [-2/100*x(1)+2/100*x(2); 2/100*x(1)-2/100*x(2)];
[t,sol] = ode45(deq1,[0 runtime],[initialx1 initialx2]);

arraysize = size(t);  % Sets time array size for the for loop.

%============ Solutions are plotted at each time step =====================

for i = 1 : max(arraysize) 
    subplot(3,1,1)
    plot(sol(iterations,1),sol(iterations,2),'.','color',[rand; rand; rand],'markersize',14,'MarkerFaceColor','b');                               
    hold on
    title(['Mixed problem Equations       t = ' num2str(t(iterations))],'fontsize',12)
    xlabel('T1','fontsize',12)
    ylabel('T2','fontsize',12)
    axis([min(sol(:,1)) max(sol(:,1)) min(sol(:,2)) max(sol(:,2))])
  
    subplot(3,1,2)
    plot(t(i,1),sol(i,1),'b.','markersize',10,'MarkerFaceColor','b')
    hold on;
    plot(t(i,1),sol(i,2),'r.','markersize',10,'MarkerFaceColor','b')
    hold on;
    axis([0 runtime 0 150])
      
    
  
    iterations = iterations + 1;   % Adds 1 to the iteration count. 
    pause(pausetime)
end

% ==== Plots time series of x or y graph depending on choice ============== 
%if choice == 1
%     subplot(3,1,2)
%     plot(t(:,1),sol(:,1),'b.','markersize',10,'MarkerFaceColor','b')
%     title(['Lotka-Volterra time series for x        run time = ' num2str(max(t)) ' seconds '],'fontsize',12)
%     xlabel('t')
%     ylabel('x')
%     axis([min(t(:,1)) max(t(:,1)) min(sol(:,1)) max(sol(:,1))])
% %elseif choice == 2
%     %title(['Lotka-Volterra time series for y        run time = ' num2str(max(t)) ' seconds '],'fontsize',12)
%     subplot(3,1,3)
%     plot(t(:,1),sol(:,2),'r.','markersize',10,'MarkerFaceColor','b')
%     title(['zoom in'],'fontsize',12)
%     xlabel('t')
%     ylabel('W')
%     axis([min(t(:,1)) max(t(:,1)) min(sol(:,2)) max(sol(:,2))])
%end

%====================== End of program ====================================