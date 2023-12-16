%% --------------- Setup --------------- %%
clc;
clear all;
close all;

%% --------------- Select Input Folder --------------- %%
PathNameInput = 'C:\Users\user\Documents\code_git\Orbit-Determination\homework-final\input-data-igs';
%PathNameInput = uigetdir(addpath(genpath(pwd)), 'Please Select Input Folder.');
FileList = dir(fullfile(PathNameInput, '*.PLT'));

%% --------------- Select Output Folder --------------- %%
PathNameOutput = 'C:\Users\user\Documents\code_git\Orbit-Determination\homework-final\output-data-igs';
%PathNameOutput = uigetdir(addpath(genpath(pwd)), 'Please Select Output Folder.');

%% --------------- Load Input Data --------------- %%
cnt_dataraw = 1;
for i = 1: numel(FileList)
    FileFid = fopen(fullfile(PathNameInput, FileList(i).name));
    disp(['User selected ', fullfile(PathNameInput, FileList(i).name)]);
    while ~feof(FileFid)
        % ---------- Split Data ---------- %
        FileTmp = (strsplit(fgetl(FileFid), ' '));
        
        % ---------- Check Raw Data ---------- %
        if size(FileTmp, 2) == 9
            index = 0;
            for j = 1: 5
                if j == 1; index = 6; end
                if j == 2; index = 3; end
                if j == 3; index = 7; end
                if j == 4; index = 8; end
                if j == 5; index = 9; end
                
                DataRaw(cnt_dataraw, j) = str2double(FileTmp(index));
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
    scatter_diff_igs(DataSat(:, 2), DataSat(:, 3: 5), sprintf('Satellite %02d', i), FigureName);
end
