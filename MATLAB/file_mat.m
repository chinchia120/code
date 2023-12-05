a = magic(3);
save my_data1.mat -ascii;
clear all;
a = load("my_data1.mat", "-ascii");

