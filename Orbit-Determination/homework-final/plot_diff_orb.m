%% --------------- Setup --------------- %%
clc;
clear all;
close all;

%% --------------- Select Input Folder --------------- %%
PathNameInput = 'C:\Users\user\Documents\code_git\Orbit-Determination\homework-final\input-data-orb';
PathNameInput = uigetdir(addpath(genpath(pwd)), 'Please Select Input Folder.');
FileList = dir(fullfile(PathNameInput, '*.PLT'));

%% --------------- Select Output Folder --------------- %%
PathNameOutput = 'C:\Users\user\Documents\code_git\Orbit-Determination\homework-final\output-data-orb';
PathNameOutput = uigetdir(addpath(genpath(pwd)), 'Please Select Output Folder.');

%% --------------- Load Input Data --------------- %%
cnt_dataraw = 1;
sat_num = 0;
for i = 1: numel(FileList)
    FileFid = fopen(fullfile(PathNameInput, FileList(i).name));
    disp(['User selected ', fullfile(PathNameInput, FileList(i).name)]);
    while ~feof(FileFid)
        % ---------- Split Data ---------- %
        FileTmp = (strsplit(fgetl(FileFid), ' '));
        
        % ---------- Check Satallite Number ---------- %
        if size(FileTmp, 2) == 6
            if FileTmp(2) == "SATELLITE" && FileTmp(3) == "NUMBER"
                sat_num = str2double(FileTmp(5));
            end
        end
        
        % ---------- Check Raw Data ---------- %
        if size(FileTmp, 2) == 7
            for j = 1: 7
                DataRaw(cnt_dataraw, 1) = sat_num;
                DataRaw(cnt_dataraw, j+1) = str2double(FileTmp(j));
            end      
        cnt_dataraw = cnt_dataraw + 1;
        end
    end
    fclose(FileFid);
end

%% --------------- Make Plot --------------- %%
for i = 1: 32
    cnt_sat = 1;
    FigureName = sprintf('%s\\Satellite%02d.png', PathNameOutput, i);
    % ---------- Select Same Satellite ---------- %
    for j = 1: size(DataRaw, 1)
        if DataRaw(j, 1) == i
            DataSat(cnt_sat, :) = DataRaw(j, :);
            cnt_sat = cnt_sat + 1;
        end
    end

    % ---------- Make Scatter ---------- %
    scatter_diff_orb(DataSat(:, 2), DataSat(:, 3: 5), sprintf('Satellite %02d', i), FigureName);
end
