clc; close all; clear; tic
% addpath ..\..\functions\
% addpath ..\..\'PCA Package Code'\
% % load("..\resultados_coagulação.mat")
% load("lambda.mat")
Data_1 = [];
Replica_1 = [];
Name_1 = [];
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
        name = categorical(string(filename(1:end-5)));
        replica=categorical(string(filename(end-4)));

        % strVector = repmat(name, size(Data_new, 2), 1);
        Replica_1=[Replica_1; replica];
        Name_1 = [Name_1; name];
        Data_1 = [Data_1; a'];
        
    end
end
% save('C&C_ON.mat',"Replica","Name","Data")

load('..\..\SISAV\SISAV.mat')

clear name
name(1:3,1)=NAME(82:84);
name(4:6,1)=NAME(79:81);
name(7:12+6,1)=NAME(67:72+6);
name=categorical(name);
clear NAME
load('..\Coagulações_triplicados\C&C_wlo.mat')

replica=[Replica;Replica_1];
x=[Data; Data_1];
NAME=[Name; name];

CeC.replica=replica;
CeC.x=x;
CeC.name=NAME;
0'