%% --------------- Setup --------------- %%
clc;
clear all;
close all;

%% --------------- Select Input brdc Folder --------------- %%
PathNameInput_brdc = 'C:\Users\user\Documents\code_git\Orbit-Determination\input-data-brdc-2010';
%PathNameInput_brdc = uigetdir(addpath(genpath(pwd)), 'Please Select Input brdc Folder.');
FileList_brdc = dir(fullfile(PathNameInput_brdc, '*.sp3'));

%% --------------- Select Input igs Folder --------------- %%
PathNameInput_igs = 'C:\Users\user\Documents\code_git\Orbit-Determination\input-data-igs-2010';
%PathNameInput_igs = uigetdir(addpath(genpath(pwd)), 'Please Select Input igs Folder.');
FileList_igs = dir(fullfile(PathNameInput_igs, '*.sp3'));

%% --------------- Select Output Folder --------------- %%
PathNameOutput = 'C:\Users\user\Documents\code_git\Orbit-Determination\output-data-trajectory-2010';
%PathNameOutput = uigetdir(addpath(genpath(pwd)), 'Please Select Output Folder.');

%% --------------- Load brdc Data --------------- %%
cnt_dataraw_brdc = 1;
for i = 1: numel(FileList_brdc)
    FileFid_brdc = fopen(fullfile(PathNameInput_brdc, FileList_brdc(i).name));
    disp(['User selected ', fullfile(PathNameInput_brdc, FileList_brdc(i).name)]);
    while ~feof(FileFid_brdc)
        FileTmp_brdc = (strsplit(fgetl(FileFid_brdc), ' '));
        
        if(size(FileTmp_brdc, 2) == 5)
            sat_num_brdc = convertStringsToChars(string(FileTmp_brdc(1)));
            if sat_num_brdc(1: 2) == 'PG'
                DataRaw_brdc{cnt_dataraw_brdc, 1} = FileTmp_brdc(1);
                DataRaw_brdc{cnt_dataraw_brdc, 2} = str2double(FileTmp_brdc(2));
                DataRaw_brdc{cnt_dataraw_brdc, 3} = str2double(FileTmp_brdc(3));
                DataRaw_brdc{cnt_dataraw_brdc, 4} = str2double(FileTmp_brdc(4));
                cnt_dataraw_brdc = cnt_dataraw_brdc + 1;
            end
        end
    end
    fclose(FileFid_brdc);
end

%% --------------- Load igs Data --------------- %%
cnt_dataraw_igs = 1;
for i = 1: numel(FileList_igs)
    FileFid_igs = fopen(fullfile(PathNameInput_igs, FileList_igs(i).name));
    disp(['User selected ', fullfile(PathNameInput_igs, FileList_igs(i).name)]);
    while ~feof(FileFid_igs)
        FileTmp_igs = (strsplit(fgetl(FileFid_igs), ' '));
          
        if(size(FileTmp_igs, 2) == 10)
            sat_num_igs = convertStringsToChars(string(FileTmp_igs(1)));
            if sat_num_igs(1: 2) == 'PG'
                DataRaw_igs{cnt_dataraw_igs, 1} = FileTmp_igs(1);
                DataRaw_igs{cnt_dataraw_igs, 2} = str2double(FileTmp_igs(2));
                DataRaw_igs{cnt_dataraw_igs, 3} = str2double(FileTmp_igs(3));
                DataRaw_igs{cnt_dataraw_igs, 4} = str2double(FileTmp_igs(4));
                cnt_dataraw_igs = cnt_dataraw_igs + 1;
            end
        end
    end
    fclose(FileFid_igs);
end

%% --------------- Make Plot --------------- %%
for i = 1: 32
    cnt_sat_brdc = 1;
    cnt_sat_igs = 1;
    FigureName = sprintf('%s\\Satellite%02d.png', PathNameOutput, i);
    TitleName = sprintf('Satellite %02d', i);
    for j = 1: max(size(DataRaw_brdc, 1), size(DataRaw_igs, 1))
        if j <= size(DataRaw_brdc, 1)
            tmp_brdc = convertStringsToChars(string(DataRaw_brdc(j, 1)));

            if str2num(tmp_brdc(3: 4)) == i
                DataSat_brdc(cnt_sat_brdc, :) = DataRaw_brdc(j, :);
                cnt_sat_brdc = cnt_sat_brdc + 1;
            end
        end

        if j <= size(DataRaw_igs, 1)
            tmp_igs = convertStringsToChars(string(DataRaw_igs(j, 1)));

            if str2num(tmp_igs(3: 4)) == i
                DataSat_igs(cnt_sat_igs, :) = DataRaw_igs(j, :);
                cnt_sat_igs = cnt_sat_igs + 1;
            end
        end
    end
    scatter_sp3(cell2mat(DataSat_brdc(:, 2: 4)), cell2mat(DataSat_igs(:, 2: 4)), TitleName, FigureName);
end
