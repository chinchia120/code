tic

A = zeros(2000, 2000);
for  i = 1: size(A, 1)
    for j = 1: size(A, 2)
        A(i, j) = i + j;
    end
end

toc

A = [1 2 3 ; ... 
    4 5 6]

