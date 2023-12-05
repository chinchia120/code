clc;
clear all;
close all;

%%% ---------- import file ---------- %%%
filename_orientation_angle = "Orientation_Angle.txt";
filename_exterior_orientation = "Exterior_Orientation.txt";

%%% ---------- check data ---------- %%%
[img, omega, phi, kappa] = textread(filename_orientation_angle ,'%s %f %f %f');
[img_, Xc, Yc, Zc, Azimuth, Elevation, Roll] = textread(filename_exterior_orientation,'%s %f %f %f %f %f %f');
rotation_matrix = zeros(size(img, 1)*3, 3);

for i = 1: size(img, 1)
    rotation_matrix(1+(i-1)*3: 3+(i-1)*3, :) = [cos(phi(i))*cos(kappa(i)), -cos(phi(i))*sin(kappa(i)), sin(phi(i));
                                                cos(omega(i))*sin(kappa(i))+sin(omega(i))*sin(phi(i))*cos(kappa(i)), cos(omega(i))*cos(kappa(i))-sin(omega(i))*sin(phi(i))*sin(kappa(i)), -sin(omega(i))*cos(phi(i));
                                                sin(omega(i))*sin(kappa(i))-cos(omega(i))*sin(phi(i))*cos(kappa(i)), sin(omega(i))*cos(kappa(i))+cos(omega(i))*sin(phi(i))*sin(kappa(i)), cos(omega(i))*cos(phi(i))];
    data(i, :) = {cell2mat(img(i)), [Xc(i), Yc(i), Zc(i)], rotation_matrix(1+(i-1)*3: 3+(i-1)*3, :)};
end

%%% ---------- match image ---------- %%%
cnt_img = 1;
for i = 1: size(img, 1)
    for j = 1 : size(img, 1)
        if i == j
            continue;
        end
        str1 = cell2mat(data(i, 1));
        str1 = str1(2: strlength(str1));
        str2 = cell2mat(data(j, 1));
        str2 = str2(2: strlength(str2));
        if strcmp(str1, str2)
            img_match(cnt_img, :) = {cell2mat(data(i, 1)), cell2mat(data(j, 1)), cell2mat(data(i, 2)), cell2mat(data(j, 2)), cell2mat(data(i, 3)), cell2mat(data(j, 3))};
            cnt_img = cnt_img + 1;
            break;
        end
    end
end
img_match = img_match(1: size(img_match, 1)/2, :);

%%% ---------- calculate ---------- %%%
for i = 1: size(img_match, 1)
    matrix = inv(cell2mat(img_match(i, 5)))*cell2mat(img_match(i, 6))
    B = inv(cell2mat(img_match(i, 5)))*transpose(cell2mat(img_match(i, 3))-cell2mat(img_match(i, 4)));
    relative_orientation(i, :) = {cell2mat(img_match(i, 1)), cell2mat(img_match(i, 2)), transpose(B), asin(matrix(1, 3)), atan(-matrix(2, 3)/matrix(3, 3)), atan(-matrix(1, 2)/matrix(1, 1))};
end

%%% ---------- output txt ---------- %%%
file = fopen('Relative_Orientation.txt', 'w');
fprintf(file, '%%%% -------------------- Relative Orientation Summary -------------------- %%%%\n\n');
fprintf(file, 'img_R\t\timg_L\t\tXc(mm)\t\t\tYc(mm)\t\t\tZc(mm)\t\t\tomega(rad)\t\t phi(rad)\t\t kappa(rad)\n');
for i = 1: size(relative_orientation, 1)
    tmp = cell2mat(relative_orientation(i, 3));
    fprintf(file, '%s\t\t%s\t\t%10.6f\t\t%10.6f\t\t%10.6f\t\t%11.8f\t\t%11.8f\t\t%11.8f\n', cell2mat(relative_orientation(i, 1)), cell2mat(relative_orientation(i, 2)), tmp(1), tmp(2), tmp(3), cell2mat(relative_orientation(i, 4)), cell2mat(relative_orientation(i, 5)), cell2mat(relative_orientation(i, 6)));
end
tmp = cell2mat(relative_orientation(:, 3));
fprintf(file, '\n\tAvg(mm/rad)\t\t\t%10.6f(mm)\t%10.6f(mm)\t%10.6f(mm)\t%11.8f(rad)%11.8f(rad)%11.8f(rad)\n', mean(tmp(:, 1)), mean(tmp(:, 2)), mean(tmp(:, 3)), mean(cell2mat(relative_orientation(:, 4))), mean(cell2mat(relative_orientation(:, 5))), mean(cell2mat(relative_orientation(:, 6))));
fprintf(file, '\tStd(mm/rad)\t\t\t%10.6f(mm)\t%10.6f(mm)\t%10.6f(mm)\t%11.8f(rad)%11.8f(rad)%11.8f(rad)\n', std(tmp(:, 1)), std(tmp(:, 2)), std(tmp(:, 3)), std(cell2mat(relative_orientation(:, 4))), std(cell2mat(relative_orientation(:, 5))), std(cell2mat(relative_orientation(:, 6))));
fprintf(file, '\n\tAvg(m/rad)\t\t\t%10.6f(m)\t%10.6f(m)\t%10.6f(m)\t%8.4f(deg)\t%7.4f(deg)\t%7.4f(deg)\n', mean(tmp(:, 1))/1000, mean(tmp(:, 2))/1000, mean(tmp(:, 3))/1000, rad2deg(mean(cell2mat(relative_orientation(:, 4)))), rad2deg(mean(cell2mat(relative_orientation(:, 5)))), rad2deg(mean(cell2mat(relative_orientation(:, 6)))));
fprintf(file, '\tStd(m/rad)\t\t\t%10.6f(m)\t%10.6f(m)\t%10.6f(m)\t%8.4f(deg)\t%7.4f(deg)\t%7.4f(deg)\n', std(tmp(:, 1))/1000, std(tmp(:, 2))/1000, std(tmp(:, 3))/1000, rad2deg(std(cell2mat(relative_orientation(:, 4)))), rad2deg(std(cell2mat(relative_orientation(:, 5)))), rad2deg(std(cell2mat(relative_orientation(:, 6)))));

%%% ---------- 95% confidence interval ---------- %%%
Xc_95_confidence_interval = sprintf('%10.6f, %10.6f',mean(tmp(:, 1)) - 2*std(tmp(:, 1)), mean(tmp(:, 1)) + 2*std(tmp(:, 1)))
Yc_95_confidence_interval = sprintf('%10.6f, %10.6f',mean(tmp(:, 2)) - 2*std(tmp(:, 2)), mean(tmp(:, 2)) + 2*std(tmp(:, 2)))
Zc_95_confidence_interval = sprintf('%10.6f, %10.6f',mean(tmp(:, 3)) - 2*std(tmp(:, 3)), mean(tmp(:, 3)) + 2*std(tmp(:, 3)))
omega_95_confidence_interval = sprintf('%11.8f, %11.8f',mean(cell2mat(relative_orientation(:, 4))) - 2*std(cell2mat(relative_orientation(:, 4))), mean(cell2mat(relative_orientation(:, 4))) + 2*std(cell2mat(relative_orientation(:, 4))))
phi_95_confidence_interval = sprintf('%11.8f, %11.8f',mean(cell2mat(relative_orientation(:, 5))) - 2*std(cell2mat(relative_orientation(:, 5))), mean(cell2mat(relative_orientation(:, 5))) + 2*std(cell2mat(relative_orientation(:, 5))))
kappa_95_confidence_interval = sprintf('%11.8f, %11.8f',mean(cell2mat(relative_orientation(:, 6))) - 2*std(cell2mat(relative_orientation(:, 6))), mean(cell2mat(relative_orientation(:, 6))) + 2*std(cell2mat(relative_orientation(:, 6))))