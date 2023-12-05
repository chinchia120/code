file = fopen('INPUT_lon_lat.txt', 'w');

for i = 119500: 25: 122500
    for j = 21500: 25: 25500
        fprintf(file, '%f\t%f\n', i/1000, j/1000);
    end
end

