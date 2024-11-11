clc; close all; clear; tic
% addpath ..\..\functions\
% addpath ..\..\'PCA Package Code'\
% % load("..\resultados_coagulação.mat")
% load("lambda.mat")
Data = [];
Replica = [];
Name = [];
% Y=[];
% Get the current working directory (the directory where the script is located)
folderPath = fullfile(pwd);
% Get a list of all folders and subfolders
folders = genpath(folderPath);
c = 1;


% Convert the semicolon-separated string into a cell array
folderList = strsplit(folders, ';');
for j = 1:numel(folderList)
    % Get a list of all XLSX files in the folder
    fileList = dir(fullfile(folderList{j}, '*.csv'));
    
    % Loop through each XLSX file and import the data
    for i = 1:numel(fileList)
        % Construct the full file path
        filePath = fullfile(folderList{j}, fileList(i).name);
        filename = fileList(i).name;
        
        % Import data from the XLSX file
        a = importfile(filename);
        
        % Extract text from the folder path
        name = categorical(string(filename(2:end-5)));
        replica=categorical(string(filename(end-4)));

        % strVector = repmat(name, size(Data_new, 2), 1);
        Replica=[Replica; replica];
        Name = [Name; name];
        Data = [Data; a'];
        
    end
end
save('C&C_wlo.mat',"Replica","Name","Data")