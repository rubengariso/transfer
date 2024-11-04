tic; clear; clc


load('data_pre_teste.mat')

addpath functions\


%% 
load("scn_phase_NIR")


idx=(find(ismember(ID, ID_teste)));
ID_TESTE=ID(idx)
x_teste=nan(length(idx),length(X_FTIR(1,:)));
y_teste=nan(length(idx),1);
for i=1:length(idx)
    x_teste(i,:)=X_FTIR(idx(i),:);
    y_teste(i)=y(idx(i));
end

idx=(find(ismember(ID, ID_tr)));
ID_TREINO=ID(idx)
x_treino=nan(length(idx),length(X_FTIR(1,:)));
y_treino=nan(length(idx),1);
for i=1:length(idx)
    x_treino(i,:)=X_FTIR(idx(i),:);
    y_treino(i)=y(idx(i));
end



save("scn_phase_FTIR", "y_treino","x_treino","ID_tr","x_teste",'y_teste','ID_teste',"ID_val","ID_cal")