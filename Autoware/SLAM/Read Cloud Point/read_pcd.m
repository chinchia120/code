%% =============== Setup =============== %%
clear all;
close all;

%% =============== Select PCD File =============== %%
pcFileName = 'Tai39_TC_193055_to_193135_@239m_1.pcd';
pcPathName = 'C:\Users\user\Documents\code_git\Autoware\direct-georeferencing\Merge Cloud Point\input-pcd\';
[pcFileName, pcPathName, ~] = uigetfile('*.pcd', 'Please Select .pcd File.');
if isequal(pcFileName, 0)
    disp('User selected Cancel');
else
    disp(['User selected ', fullfile(pcPathName, pcFileName)]);
end

if  pcFileName == 0 
    disp('No file is selected.'); 
    return; 
end

file = [pcPathName, pcFileName];
ptCloud = pcread(file);

%% =============== Show PCD =============== %%
figure;pcshow(ptCloud);
