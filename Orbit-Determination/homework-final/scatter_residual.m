%% --------------- Function of Residual Scatter --------------- %%
function output = scatter_residual(X_, Y_, title_, figurename_)
    % ---------- Subplot 1 ---------- %
    subplot(3, 1, 1);
    scatter(X_, Y_(:, 1), '.', 'b');

    title([title_, ' of Direction X'], 'FontSize', 16);
    xlabel('Time', 'FontSize', 14);
    ylabel('Residual (m)', 'FontSize', 14);
    
    % ---------- Subplot 2 ---------- %
    subplot(3, 1, 2);
    scatter(X_, Y_(:, 2), '.', 'r');

    title([title_, ' of Direction Y'], 'FontSize', 16);
    xlabel('Time', 'FontSize', 14);
    ylabel('Residual (m)', 'FontSize', 14);
    
    % ---------- Subplot 3 ---------- %
    subplot(3, 1, 3);
    scatter(X_, Y_(:, 3), '.', 'k');

    title([title_, ' of Direction Z'], 'FontSize', 16);
    xlabel('Time', 'FontSize', 14);
    ylabel('Residual (m)', 'FontSize', 14);
    
    % ---------- Save Figure ---------- %
    set(gcf, 'position', [0, 0, 1250, 625]);
    saveas(gcf, figurename_);
    output = [];
end
