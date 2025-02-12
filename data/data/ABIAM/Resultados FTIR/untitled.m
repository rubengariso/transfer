clc; close all; clear; tic
addpath ..\..\functions\
addpath ..\..\'PCA Package Code'\
load("..\resultados_coagulação.mat")
load("lambda.mat")
Data = [];
Name = [];
Sample = [];
Y=[];
% Get the current working directory (the directory where the script is located)
folderPath = fullfile(pwd);
% Get a list of all folders and subfolders
folders = genpath(folderPath);
c = 1;

coag = data.coagula;
ID = categorical(data.ID);  

% Convert the semicolon-separated string into a cell array
folderList = strsplit(folders, ';');
for j = 1:numel(folderList)
    % Get a list of all XLSX files in the folder
    fileList = dir(fullfile(folderList{j}, '*.xlsx'));
    
    % Loop through each XLSX file and import the data
    for i = 1:numel(fileList)
        % Construct the full file path
        filePath = fullfile(folderList{j}, fileList(i).name);
        filename = fileList(i).name;
        
        % Import data from the XLSX file
        [Data_new, Name_new] = importfile(filePath);
        
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



incremento=0.1

x=2- log10(Data + incremento);

Sample_string = string(Sample);
ABIAN.spectra=Data;
ABIAN.sample=Sample;
ABIAN.lambda=lambeda;
% Define the color that should be used
rnd = colormap('lines');
symbo = {'*', 's', 'o', '<', 'x', 'v'};
close


Y_unique=unique(Y);
figure
plot(lambeda,Data)
% xlabel(sprintf('Score #1 (%.1f%%)', CP1));
% ylabel(sprintf('Score #2 (%.1f%%)', CP2));
set(gca, 'FontSize', 16)
figure
hold on
for i=1:5
    idx=ismember(Y,Y_unique(i));
    plot(lambeda,(x(idx,:)) );
end
legend(Sample(idx), 'Location', 'BestOutside')
set(gca, 'FontSize', 16)




figure
hold on
for i=5
    idx=ismember(Y,Y_unique(i));
    plot(lambeda,(x(idx,:)) );
end
legend(Sample(idx), 'Location', 'BestOutside')
set(gca, 'FontSize', 16)

% ind_1=find(Y==1 & x(:,314)<2,1,'first');
ind_0=find(ismember(Sample,Sample(end-3*1+1)),1,'first');
ind_1=find(ismember(Sample,Sample(end-3*2+1)),1,'first');
idx=[ind_0 ind_1];
color_plot=[0 ,176/255, 113/255; 1,0,0];
for i=1:2;
    figure
% idx=ismember(Sample,Sample(end-3*i+1));
% plot(lambeda,mean(x(idx,:)),'LineWidth',8,'Color',color_plot(i,:))
% [0 ,176/255, 113/255]);
plot(lambeda,x(idx(i),:),'LineWidth',5,'Color','k')
ylim([0 3])
set(gca,'XColor', 'none','YColor','none')
end