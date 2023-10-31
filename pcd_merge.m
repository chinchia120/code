%% ========== Setup ========== %%
clear all;
close all;

%% ========== Select First PCD File ========== %%
%[pcFileName1, pcPathName1, ~] = uigetfile('*.pcd', 'Please Select .pcd File.');
pcFileName1 = 'Tai39_TC_193055_to_193135_@239m_1.pcd';
pcPathName1 = 'C:\Users\user\Documents\code_git\Autoware\direct-georeferencing\Merge Cloud Point\';
if isequal(pcFileName1, 0)
    disp('User selected Cancel');
else
    disp(['User selected ', fullfile(pcPathName1, pcFileName1)]);
end

if  pcFileName1 == 0 
    disp('No file is selected.'); 
    return; 
end

file1 = [pcPathName1, pcFileName1];
ptCloud1 = pcread(file1);

%% ========== Select Second PCD File ========== %%
%[pcFileName2, pcPathName2, ~] = uigetfile('*.pcd', 'Please Select .pcd File.');
pcFileName2 = 'Tai39_TC_193122_to_193199_@724m_1.pcd';
pcPathName2 = 'C:\Users\user\Documents\code_git\Autoware\direct-georeferencing\Merge Cloud Point\';
if isequal(pcFileName2, 0)
    disp('User selected Cancel');
else
    disp(['User selected ', fullfile(pcPathName2, pcFileName2)]);
end

if  pcFileName2 == 0 
    disp('No file is selected.'); 
    return; 
end

file2 = [pcPathName2, pcFileName2];
ptCloud2 = pcread(file2);

%% ========== Select Output PCD Folder ========== %%
%outputPathName = uigetdir(matlabroot, 'Please Select Output Folder.');
outputPathName = 'C:\Users\user\Documents\code_git\Autoware\direct-georeferencing\Merge Cloud Point\Outoput-PCD';
outputFileName = 'merge.pcd';
file_merge = [outputPathName, '\', outputFileName];

%% ========== Merge PCD ========== %%
ptCloud_merge = pcmerge(ptCloud1, ptCloud2, 0.0001);

%% ========== Show PCD ========== %%
figure;pcshow(ptCloud1);
figure;pcshow(ptCloud2);
figure;pcshow(ptCloud_merge);

%% ========== Output PCD ========== %%
pcwrite(ptCloud_merge, file_merge);
