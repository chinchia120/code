A = 20;
B = int8(20);
whos;

str = 'abcdea';
str(3);
'a' == str;
str == 'a';
str(str == 'a') = 'Z';

j = strlength(str) + 1;
for i = 1: strlength(str)
    s(i) = str(j - i);
end
s
