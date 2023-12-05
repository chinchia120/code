%% =============== Setup =============== %%
clc;
clear all;

%% =============== Select Folder1 =============== %%
PathName1 = 'C:\Users\user\Documents\code_git\TA\Photogrammetry\data-output-student\F64101074';
PathName1 = uigetdir(addpath(genpath(pwd)), 'Please Select First Folder.');
if isequal(PathName1, 0)
    disp('User selected Cancel');
    return;
else
    disp(['User selected ', PathName1]);
end
FileList1 = dir(fullfile(PathName1, '*.out'));

%% =============== Select Folder2 =============== %%
PathName2 = 'C:\Users\user\Documents\code_git\TA\Photogrammetry\data-output-check\F64101074';
PathName2 = uigetdir(addpath(genpath(pwd)), 'Please Select Second Folder.');
if isequal(PathName2, 0)
    disp('User selected Cancel');
    return;
else
    disp(['User selected ', PathName2]);
end
FileList2 = dir(fullfile(PathName2, '*.out'));
disp([sprintf('<< =============== Check Folder =============== >>\n')])
%% =============== Compare Process =============== %%
for i = 1: numel(FileList1)
    FileFid1 = fopen(fullfile(PathName1, FileList1(i).name));
    disp(['User selected ', fullfile(PathName1, FileList1(i).name)]);
    
    % ===== Loading File ===== %
    cnt1 = 1;
    cnt2 = 1;
    for j = 1: numel(FileList2)
        if(strcmp(FileList1(i).name, FileList2(j).name))
            FileFid2 = fopen(fullfile(PathName2, FileList2(j).name));
            disp(['User selected ', fullfile(PathName2, FileList2(j).name)]);
            
            while ~feof(FileFid1)
            Filetmp1 = fgetl(FileFid1);
            DataRaw1(cnt1, :) = {Filetmp1};
            cnt1 = cnt1 + 1;
            end

            while ~feof(FileFid2)
            Filetmp2 = fgetl(FileFid2);
            DataRaw2(cnt2, :) = {Filetmp2};
            cnt2 = cnt2 + 1;
            end

            break;
        end
    end

    % ===== Compare File ===== %
    flag = 1;
    for j = 1: min(size(DataRaw1, 1), size(DataRaw2, 1))
        tmp1 = cell2mat(DataRaw1(j, :));
        tmp2 = cell2mat(DataRaw2(j, :));
        if ~strcmp(tmp1, tmp2)
            disp(tmp1);
            disp(tmp2);
            flag = 0;
        end
    end

    if flag == 0; disp([sprintf('<< ===============  ERROR  =============== >>\n')]); end
    if flag == 1; disp([sprintf('<< =============== CORRECT =============== >>\n')]); end

    % ===== Close File ===== %
    fclose(FileFid1);
    fclose(FileFid2);
    clearvars DataRaw1 DataRaw2;
end
