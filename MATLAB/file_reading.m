
fid = fopen("data.txt", "r");
i = 1;
while ~feof(fid)
    name(i, :) = fscanf(fid, '%c', 1);
    year(i) = fscanf(fid, '%d', 1);
    month(i) = fscanf(fid, '%d', 1);
    day(i) = fscanf(fid, '%d\n', 1);
    i = i + 1;
end
fclose(fid);
data = {name ; year; month; day}