clc; close all; clear; tic
% addpath ..\..\functions\
% addpath ..\..\'PCA Package Code'\
% % load("..\resultados_coagulação.mat")
% load("lambda.mat")
Data = [];
Name = [];
Sample = [];
Y=[];
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
        extractedText = categorical(string(folderList{j}(length(folderPath) + 2:end)));
        idx=ismember(ID, extractedText);
        if any(idx~=0)
            y=coag (idx);
        switch y
            case 'Coagula'
                Y_text=1;
            case 'Não coagula' 
                Y_text=0;
            case 'Coagula - ambíguo' 
                Y_text=2;
            case 'Não coagula - ambíguo'
                Y_text=3;
            otherwise
                puta
                end

        else
            Y_text=5;
        end
            
        % Append new data and sample information
        strVector = repmat(extractedText, size(Data_new, 2), 1);
        Y=[Y; repmat(Y_text, size(Data_new, 2), 1)];
        Sample = [Sample; strVector];
        Data = [Data; Data_new'];
        Name = [Name; Name_new(1, :)];
    end
end
