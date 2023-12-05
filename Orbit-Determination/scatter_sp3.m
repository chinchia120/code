%% --------------- Function of sp3 Scatter --------------- %%
function output = scatter_sp3(data1_, data2_, title_, figurename_)
    % ---------- Setup ---------- %
    [X, Y, Z] = sphere;

    % ---------- Subplot 1 ---------- %
    subplot(1, 3, 1);

    scatter3(data1_(:, 1), data1_(:, 2), data1_(:, 3), '+', 'b');
    hold on;
    surf(X*6371, Y*6371, Z*6371);
    hold on;

    title([title_, ' of brdc Data'], 'FontSize', 16);
    xlabel('X - axis', 'FontSize', 14);
    ylabel('Y - axis', 'FontSize', 14);
    zlabel('Z - axis', 'FontSize', 14);
    hold off;

    % ---------- Subplot 2 ---------- %
    subplot(1, 3, 2);

    scatter3(data2_(:, 1), data2_(:, 2), data2_(:, 3), 'o', 'r');
    hold on;
    surf(X*6371, Y*6371, Z*6371);
    hold on;

    title([title_, ' of igs Data'], 'FontSize', 16);
    xlabel('X - axis', 'FontSize', 14);
    ylabel('Y - axis', 'FontSize', 14);
    zlabel('Z - axis', 'FontSize', 14);
    hold off;

    % ---------- Subplot 3 ---------- %
    subplot(1, 3, 3);

    scatter3(data1_(:, 1), data1_(:, 2), data1_(:, 3), '+', 'b');
    hold on;
    scatter3(data2_(:, 1), data2_(:, 2), data2_(:, 3), 'o', 'r');
    hold on;
    surf(X*6371, Y*6371, Z*6371);
    hold on;

    title([title_, ' of Merge Data'], 'FontSize', 16);
    xlabel('X - axis', 'FontSize', 14);
    ylabel('Y - axis', 'FontSize', 14);
    zlabel('Z - axis', 'FontSize', 14);
    hold off;
    
    % ---------- Save Figure ---------- %
    set(gcf, 'position', [0, 0, 1250, 625]);
    saveas(gcf, figurename_);
    output = [];
end

