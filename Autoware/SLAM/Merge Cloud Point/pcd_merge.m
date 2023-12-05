%% =============== Setup =============== %%
clc;
clear all;

%% =============== Select First PCD File =============== %%
pcFileName1 = 'vlp16_DG_06291_09740_.pcd';
pcPathName1 = 'C:\Users\user\Documents\code_git\Autoware\SLAM\Merge Cloud Point\input-pcd\';
[pcFileName1, pcPathName1, ~] = uigetfile('*.pcd', 'Please Select .pcd File.');
if isequal(pcFileName1, 0)
    disp('User selected Cancel');
    return;
else
    disp(['User selected ', fullfile(pcPathName1, pcFileName1)]);
end

file1 = [pcPathName1, pcFileName1];
ptCloud1 = pcread(file1);

%% =============== Select Second PCD File =============== %%
pcFileName2 = 'vlp16_DG_12870_14021_.pcd';
pcPathName2 = 'C:\Users\user\Documents\code_git\Autoware\SLAM\Merge Cloud Point\input-pcd\';
[pcFileName2, pcPathName2, ~] = uigetfile('*.pcd', 'Please Select .pcd File.');
if isequal(pcFileName2, 0)
    disp('User selected Cancel');
    return;
else
    disp(['User selected ', fullfile(pcPathName2, pcFileName2)]);
end

file2 = [pcPathName2, pcFileName2];
ptCloud2 = pcread(file2);

%% =============== Select Output PCD Folder =============== %%
outputPathName = 'C:\Users\user\Documents\code_git\Autoware\SLAM\Merge Cloud Point\output-pcd';
%outputPathName = uigetdir(addpath(genpath(pwd)), 'Please Select Output Folder.');

File1_spt = strsplit(pcFileName1, '_');
File2_spt = strsplit(pcFileName2, '_');
File1_time = [str2num(cell2mat(File1_spt(3))), str2num(cell2mat(File1_spt(4)))];
File2_time = [str2num(cell2mat(File2_spt(3))), str2num(cell2mat(File2_spt(4)))];
outputFileName = sprintf('%s_%s_%05d_%05d_.pcd', cell2mat(File1_spt(1)), cell2mat(File1_spt(2)), min(File1_time(1), File2_time(1)), max(File1_time(2), File2_time(2)));
file_merge = [outputPathName, '\', outputFileName];

%% =============== Merge PCD =============== %%
ptCloud_merge = pcmerge(ptCloud1, ptCloud2, 0.0001);

%% =============== Output PCD =============== %%
pcwrite(ptCloud_merge, file_merge);

%% =============== Show PCD =============== %%
%figure;pcshow(ptCloud1);
%figure;pcshow(ptCloud2);
%figure;pcshow(ptCloud_merge);
