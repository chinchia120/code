%% --------------- Setup --------------- %%
clc;
clear all;
close all;

%% --------------- Select Input Folder --------------- %%
PathNameInput = 'C:\Users\user\Documents\code_git\Orbit-Determination\input-data-residual-2010';
PathNameInput = uigetdir(addpath(genpath(pwd)), 'Please Select Input Folder.');
FileList = dir(fullfile(PathNameInput, '*.txt'));

%% --------------- Select Output Folder --------------- %%
PathNameOutput = 'C:\Users\user\Documents\code_git\Orbit-Determination\output-data-residual-2010';
PathNameOutput = uigetdir(addpath(genpath(pwd)), 'Please Select Output Folder.');

%% --------------- Load Input Data --------------- %%
cnt_dataraw = 1;
for i = 1: numel(FileList)
    FileFid = fopen(fullfile(PathNameInput, FileList(i).name));
    disp(['User selected ', fullfile(PathNameInput, FileList(i).name)]);
    while ~feof(FileFid)
        FileTmp = (strsplit(fgetl(FileFid), ' '));
        if(size(FileTmp, 2) == 14)
            index = 0;
            for j = 1: 5
                if j == 1; index = 03; end
                if j == 2; index = 05; end
                if j == 3; index = 10; end
                if j == 4; index = 11; end
                if j == 5; index = 12; end
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
    for j = 1: size(DataRaw, 1)
        if DataRaw(j, 1) == i
            DataSat(cnt_sat, :) = DataRaw(j, :);
            cnt_sat = cnt_sat + 1;
        end
    end
    scatter_residual(DataSat(:, 2), DataSat(:, 3: 5), sprintf('Satellite %02d', i), FigureName);
end
