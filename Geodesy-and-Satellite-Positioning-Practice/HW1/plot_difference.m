close all;
clear all;

show_plot('EGM96-EGM2008.txt', 'EGM96-EGM2008', 'EGM96-EGM2008.png');
show_plot('EGM96-Twgeoid.txt', 'EGM96-Twgeoid', 'EGM96-Twgeoid.png');
show_plot('EGM96-Twgeoid_hybrid.txt', 'EGM96-Twgeoid hybrid', 'EGM96-Twgeoid_hybrid.png');
show_plot('EGM2008-Twgeoid.txt', 'EGM2008-Twgeoid', 'EGM2008-Twgeoid.png');
show_plot('EGM2008-Twgeoid_hybrid.txt', 'EGM2008-Twgeoid hybrid', 'EGM2008-Twgeoid_hybrid.png');
show_plot('Twgeoid-Twgeoid_hybrid.txt', 'Twgeoid-Twgeoid hybrid', 'Twgeoid-Twgeoid_hybrid.png');

function output = show_plot(data_name, title_name, img_name)
    x = load(data_name);
    [m, n] = size(x);
    nx = 1;
    ny = 0;
    limx = (122.5-119.5)/0.025 + 1;
    limy = (25.5-21.5)/0.025 + 1;
    for i = 1: m
        ny = ny + 1;
        im(limy - ny + 1, nx) = x(i, 3);
        if(ny == limy)
            ny = 0;
            nx = nx + 1;
        end
    end

    imagesc(im);
    title(title_name);
    set(gca, 'YTick', [1, (1+limy) / 2, limy], 'YTicklabel', {'25.75', '23.50', '21.25'});
    set(gca, 'XTick', [1, (1+limx) / 2, limx], 'XTicklabel', {'119.25', '121.00', '122.75'});
    colorbar;
    saveas(gcf, img_name);

    output = [];
end

