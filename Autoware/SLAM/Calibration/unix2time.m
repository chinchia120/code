function t = unix2time(unix_t)
%UNIX2TIME Summary of this function goes here
%   Detailed explanation goes here
    t.time = fix(unix_t);
    t.sec = unix_t - fix(unix_t);
    return;
end

