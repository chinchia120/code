file = fopen('INPUT_lat_lon.txt', 'w');

for i = 119500: 25: 122500
    for j = 21500: 25: 25500
        fprintf(file, '%f\t%f\n', j/1000, i/1000);
    end
end

