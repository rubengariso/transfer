function [ID_teste, ID_tr,ID_cal,ID_val] = defin_grup(ID0,ID1,frac, C_max)
%DEFIN_GRUP Summary of this function goes here
%   Detailed explanation goes here
% ID_teste -as amostras que sao de teste 
ID_1=unique(ID1);
ID_0=unique(ID0);

%C=100; % numero de conjuntos calibração / validação
%% conjunto de teste (20%)
x_max=min(length(ID_0),length(ID_1)); % para garantir que as amostras estao balanciadas 

N_test = round(frac* x_max);


ID_0_teste= ID_0(randperm (length(ID_0), N_test));
ID_1_teste= ID_1(randperm (length(ID_1), N_test));
ID_teste= union(ID_0_teste, ID_1_teste);


%% conjunto de treino 

id_0_tr=ID_0(~ismember(ID_0, ID_0_teste)); 
id_1_tr=ID_1(~ismember(ID_1, ID_1_teste));
ID_tr= union(id_0_tr, id_1_tr);


x_max=min(length(id_0_tr),length(id_1_tr));

N_calibr=round(0.8 * x_max);
N_val= x_max - N_calibr;

for c=1:C_max
    ID_0_cal=id_0_tr(randperm(length(id_0_tr), N_calibr));
    ID_0_v=id_0_tr(~ismember(id_0_tr,ID_0_cal));
    ID_0_v=ID_0_v(randperm(length(ID_0_v), N_val));
    ID_1_cal=id_1_tr(randperm(length(id_1_tr), N_calibr ));
    ID_1_v=id_1_tr(~ismember(id_1_tr,ID_1_cal));
    ID_1_v=ID_1_v(randperm(length(ID_1_v), N_val)); 
    id_val=union(ID_0_v,ID_1_v);
   id_cal= union(ID_0_cal, ID_1_cal);

    ID_val(c,:)=(find(ismember(ID_tr, id_val)));
   ID_cal(c,:)=(find(ismember(ID_tr, id_cal)));

end
end

