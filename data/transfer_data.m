clc; close all; clear; tic

load('Coagulacoes_FTIR_C&C_24\Amostras Sogilub_triplicados\C&C.mat')

load('ABIAM\Resultados FTIR\ABIAN.mat')

load('Carmona\carmona.mat')

load("SISAV\SISAV.mat")


addpath ..\matalab_KEYS_funcion\'PCA Package Code'
addpath ..\matalab_KEYS_funcion\functions\
%% conver dados por causa dos erros de formatos (aka experimentalistas sao m*rda)

load("morte_aos_experimentalistas.mat")
SISAV.spectra=x;
SISAV.name=categorical(NAME);
SISAV.lambda=lambda*10^(-6);
CeC.lambda=CeC.lambda*10^(-2);
axis_transfer=tranfer_v2;


SISAV.spectra=100*SISAV.spectra;
CeC.x=10^-2*CeC.x;



for i=1:size(axis_transfer,1)
    idx_A=ismember(ABIAN.sample, axis_transfer(i,1));
    ABIAN.sample(idx_A)=axis_transfer(i,2);
    idx_A=ismember(CARMONA.sample, axis_transfer(i,1));
    CARMONA.sample(idx_A)=axis_transfer(i,2);
    idx_A=ismember(CeC.name, axis_transfer(i,1));
    CeC.name(idx_A)=axis_transfer(i,2);
    idx_A=ismember(SISAV.name, axis_transfer(i,1));
    SISAV.name(idx_A)=axis_transfer(i,2);
end



%%  definir a gama de lambdas a usar;

lam_alta=[CARMONA.lambda(1);SISAV.lambda(1);ABIAN.lambda(1);CeC.lambda(1)];          lamba_min=min(lam_alta);
lam_baixa=[CARMONA.lambda(end);SISAV.lambda(end);ABIAN.lambda(end);CeC.lambda(end)]; lamba_max=max(lam_baixa);

index_CeC=find(CeC.lambda<=lamba_min & CeC.lambda>=lamba_max);
lambda_ref=CeC.lambda(index_CeC)';

%% espectros estimados

x_carmona   = spline(CARMONA.lambda,CARMONA.spectrum,lambda_ref);
x_SISAV     = spline(SISAV.lambda,SISAV.spectra,lambda_ref);
x_CeC       = spline(CeC.lambda,CeC.x,lambda_ref);
x_ABIAN     = spline(ABIAN.lambda,ABIAN.spectra,lambda_ref);


%% reorganizar os espectro de modo a ter as amostras totais


ID_global=[SISAV.name; ABIAN.sample; CARMONA.sample; CeC.name];
ID_unique=unique(ID_global);

idx_A=ismember(ID_unique,ABIAN.sample); idx_B=ismember(ID_unique,CARMONA.sample);
idx_C=ismember(ID_unique,CeC.name);     idx_D=ismember(ID_unique,SISAV.name);

idx=[idx_A, idx_B, idx_C, idx_D];

ID_comum=ID_unique(sum(idx,2)==4);
% ID_comum=sum(idx,2)==4;
% ID_not_comum=ID_unique(~is_comum)

x_ABIAN_PCA=zeros(3*size(ID_comum,1),size(lambda_ref,2));
x_CARMONA_PCA=zeros(3*size(ID_comum,1),size(lambda_ref,2));
x_CeC_PCA=zeros(3*size(ID_comum,1),size(lambda_ref,2));
x_SISAV_PCA=zeros(3*size(ID_comum,1),size(lambda_ref,2));
x_ABIAN_transfer=zeros(3*size(ID_comum,1),size(ABIAN.spectra,2));
x_CARMONA_transfer=zeros(3*size(ID_comum,1),size(CARMONA.spectrum,2));
x_CeC_transfer=zeros(3*size(ID_comum,1),size(CeC.x,2));
x_SISAV_transfer=zeros(3*size(ID_comum,1),size(SISAV.spectra,2));
IDX_start=1:3:size(x_SISAV_PCA,1);
IDX_END=3:3:size(x_SISAV_PCA,1);

for i=1:size(ID_comum,1)
    idx_A=ismember(ABIAN.sample, ID_comum(i));
    x_ABIAN_PCA(IDX_start(i):IDX_END(i),:)=x_ABIAN(idx_A,:);
    x_ABIAN_transfer(IDX_start(i):IDX_END(i),:)=ABIAN.spectra(idx_A,:);
    idx_A=ismember(CARMONA.sample, ID_comum(i));
    x_CARMONA_PCA(IDX_start(i):IDX_END(i),:)=x_carmona(idx_A,:);
    x_CARMONA_transfer(IDX_start(i):IDX_END(i),:)=CARMONA.spectrum(idx_A,:);
    idx_A=ismember(CeC.name, ID_comum(i));
    x_prov=x_CeC(idx_A,:);
    x_CeC_PCA(IDX_start(i):IDX_END(i),:)=x_prov(1:3,:);
    x_prov=CeC.x(idx_A,:);
    x_CeC_transfer(IDX_start(i):IDX_END(i),:)=x_prov(1:3,:);
    idx_A=ismember(SISAV.name, ID_comum(i));
    x_SISAV_PCA(IDX_start(i):IDX_END(i),:)=x_SISAV(idx_A,:);
    x_SISAV_transfer(IDX_start(i):IDX_END(i),:)=SISAV.spectra(idx_A,:);
    ID_final_dos_finais(IDX_start(i):IDX_END(i),:)=ID_comum(i);
end

figure
subplot(2,2,1)
plot(lambda_ref,x_ABIAN_PCA)
xlabel('\lambda (cm^{-1})')
ylabel('T(%)')
title('ABIAN')
set(gca, 'FontSize', 16)
axis([500 4000 0 90])
subplot(2,2,2)
plot(lambda_ref,x_CARMONA_PCA)
xlabel('\lambda (cm^{-1})')
ylabel('T(%)')
title('CARMONA')
axis([500 4000 0 90])
set(gca, 'FontSize', 16)
subplot(2,2,3)
plot(lambda_ref,x_CeC_PCA)
xlabel('\lambda (cm^{-1})')
ylabel('T(%)')
title('C&C')
axis([500 4000 0 90])
set(gca, 'FontSize', 16)
subplot(2,2,4)
plot(lambda_ref,x_SISAV_PCA)
xlabel('\lambda (cm^{-1})')
ylabel('T(%)')
title('SISAV')
set(gca, 'FontSize', 16)
axis([500 4000 0 90])

figure
subplot(2,2,1)
plot(ABIAN.lambda,x_ABIAN_transfer)
xlabel('\lambda')
ylabel('T')
title('ABIAN')
set(gca, 'FontSize', 16)
axis([500 4000 0 90])
subplot(2,2,2)
plot(CARMONA.lambda,x_CARMONA_transfer)
xlabel('\lambda')
ylabel('T')
title('CARMONA')
set(gca, 'FontSize', 16)
axis([500 4000 0 90])
subplot(2,2,3)
plot(CeC.lambda,x_CeC_transfer)
xlabel('\lambda')
ylabel('T')
title('C & C')
set(gca, 'FontSize', 16)
axis([500 4000 0 90])
subplot(2,2,4)
plot(SISAV.lambda,x_SISAV_transfer)
xlabel('\lambda')
ylabel('T')
title('SISAV')
set(gca, 'FontSize', 16)
axis([500 4000 0 90])



x_ABIAN=x_ABIAN_PCA -mean(x_ABIAN_PCA);
x_CARMONA=x_CARMONA_PCA- mean(x_CARMONA_PCA);
x_CeC=x_CeC_PCA-mean(x_CeC_PCA);
x_SISAV=x_SISAV_PCA-mean(x_SISAV_PCA);

x_pca=[x_ABIAN;x_CARMONA;x_CeC;x_SISAV];
alpha=1/100;
kSelect='kSelectCpv';
[ model_PCA ] = pcaModel(x_pca, kSelect, 'datamat' );



% Determe PC of each 
T_ABIAN_1   =   x_ABIAN     * model_PCA.P(:,1:model_PCA.k);
T_CARMONA_1 =   x_CARMONA   * model_PCA.P(:,1:model_PCA.k);
T_CeC_1     =   x_CeC       * model_PCA.P(:,1:model_PCA.k);
T_SISAV_1   =   x_SISAV     * model_PCA.P(:,1:model_PCA.k);


figure
hold on
plot(T_ABIAN_1(:,1),T_ABIAN_1(:,2),'b.')
plot(T_CARMONA_1(:,1),T_CARMONA_1(:,2),'r*')
plot(T_CeC_1(:,1),T_CeC_1(:,2),'ks')
plot(T_SISAV_1(:,1),T_SISAV_1(:,2),'go')
hold off
xlabel('PC1')
ylabel('PC2')
legend('ABIAN','CARMONA','CeC','SISAV')
set(gca, 'FontSize', 16)
axis([-2000 1000 -1000  500 ])
box on
% model_PCA.ucl.sT=limitT(alpha/2,model_PCA.k);
% model_PCA.ucl.sE=limitQ(model_PCA.L,model_PCA.k,alpha/2);
% [ stat_pca_A ] = pcaMonitor(x_ABIAN, model_PCA );
% [ stat_pca_CA ] = pcaMonitor(x_CARMONA, model_PCA );
% [ stat_pca_CC ] = pcaMonitor(x_CeC, model_PCA );
% [ stat_pca_S ] = pcaMonitor(x_SISAV, model_PCA );

%% direct standardization DS
load('training_transfer.mat')
% X_Ref=X_i*F


% x_CeC_transfer=x_CeC_transfer-mean(x_CeC_transfer);
% x_ABIAN_transfer=x_ABIAN_transfer-mean(x_ABIAN_transfer);
% x_CARMONA_transfer=x_CARMONA_transfer-mean(x_CARMONA_transfer);
% x_SISAV_transfer=x_SISAV_transfer-mean(x_SISAV_transfer);
% X_REF is the domain we want to project to
%X_new is the original domain
% [F_DS_ABIAN,model_ABIAN,RMSE_ABIAN,Yhat_ABIAN]          = direct_standardization(x_CeC_transfer,x_ABIAN_transfer);
% [F_DS_CARMONA,model_CARMONA,RMSE_CARMONA,Yhat_CARMONA]  = direct_standardization(x_CeC_transfer,x_CARMONA_transfer);
% [F_DS_SISAV,model_SISAV,RMSE_SISAV,Yhat_SISAV]          = direct_standardization(x_CeC_transfer,x_SISAV_transfer);
% 
% 
% 
RMSE_DS=[RMSE_ABIAN, RMSE_CARMONA, RMSE_SISAV];


name=["ABIAN-DS","CARMONA-DS","SISAV-DS"];
figure
boxplot(RMSE_DS,name)
ylabel('RMSE')
set(gca, 'FontSize', 16)
% 
% 
% 
T_ABIAN_R_DS   =   Yhat_ABIAN(:,index_CeC)     * model_PCA.P(:,1:model_PCA.k);
T_CARMONA_R_DS =   Yhat_CARMONA(:,index_CeC)   * model_PCA.P(:,1:model_PCA.k);
% T_CeC_R_DS     =   x_CeC_transfer(:,index_CeC) * model_PCA.P(:,1:model_PCA.k);
T_SISAV_R_DS   =   Yhat_SISAV(:,index_CeC)     * model_PCA.P(:,1:model_PCA.k);

figure
hold on
plot(T_ABIAN_R_DS(:,1),T_ABIAN_R_DS(:,2),'b.')
plot(T_CARMONA_R_DS(:,1),T_CARMONA_R_DS(:,2),'r*')
plot(T_CeC_1(:,1),T_CeC_1(:,2),'ks')
plot(T_SISAV_R_DS(:,1),T_SISAV_R_DS(:,2),'go')
hold off
xlabel('PC1')
ylabel('PC2')
legend('ABIAN','CARMONA','CeC','SISAV','Location','southeast')
set(gca, 'FontSize', 16)
axis([-2000 1000 -1000  500 ])
box on
% 

toc
%% piecewise direct standardization

% [F_ABIAN_PDS,RMSE_ABIAN_PDS,Yhat_ABIAN_PDS] = piecewise_direct_standardization(x_ABIAN_transfer,x_CeC_transfer,10);
% [F_CARMONA_PDS,RMSE_CARMONA_PDS,Yhat_CARMONA_PDS] = piecewise_direct_standardization(x_CARMONA_transfer,x_CeC_transfer,10);
% [F_SISAV_PDS,RMSE_SISAV_PDS,Yhat_SISAV_PDS] = piecewise_direct_standardization(x_SISAV_transfer,x_CeC_transfer,100);
% 
% RMSE=[RMSE_SISAV_PDS,RMSE_CARMONA_PDS,RMSE_ABIAN_PDS];
% name=["SISAV-PDS","CARMONA-PDS","ABIAN-PDS"];
% figure
% boxplot(RMSE,name)
% ylabel('RMSE')
% set(gca, 'FontSize', 16)
% i_max=size(lambda_ref,2);
% windoh=10;
% puta
% 
% 
% for i=1:i_max
%     j_min=i-windoh;
%     j_max=i+windoh;
%     if j_min<1
%         j_min=1;
%     end
%     if j_max>i_max
%         j_max=i_max;
%     end
%     [ ~,~,~,~,BETA, ~ ] = plsModel(x_ABIAN(:,j_min:j_max)   , x_CeC(:,i) ,size(j_min:j_max,2) );
%     F_ABIAN_PSD(j_min:j_max,i)=BETA.Y;
%     % size(x_ABIAN(:,j_min:j_max))
%     % size(F_ABIAN_PSD(j_min:j_max,i))
%     % x_reformed=x_ABIAN(:,j_min:j_max)*F_ABIAN_PSD(j_min:j_max,i);
%     [ ~,~,~,~,BETA, ~ ] = plsModel(x_CARMONA(:,j_min:j_max) , x_CeC(:,i) ,size(j_min:j_max,2) );
%     F_CARMONA_PSD(j_min:j_max,i)=BETA.Y;
%     [ ~,~,~,~,BETA, ~ ] = plsModel(x_SISAV(:,j_min:j_max)   , x_CeC(:,i) ,size(j_min:j_max,2) );
%     F_SISAV_PSD(j_min:j_max,i)=BETA.Y;
% end
% % (i,j_min:j_max)(j_min:j_max,i)(i,j_min:j_max)
% x_hat_abian=x_ABIAN*F_ABIAN_PSD;
% x_hat_car=x_CARMONA*F_CARMONA_PSD;
% x_hat_sis=x_SISAV*F_SISAV_PSD;
% 
% RMSE_sis=sqrt(sum( (x_CeC-x_hat_sis).^2,2 ));
% RMSE_ab=sqrt(sum( (x_CeC-x_hat_abian).^2,2 ));
% RMSE_car=sqrt(sum( (x_CeC-x_hat_car).^2,2 ));





RMSE_PDS=[RMSE_ABIAN_PDS, RMSE_CARMONA_PDS, RMSE_SISAV_PDS];
name=["ABIAN-PDS","CARMONA-PDS","SISAV-PDS"];
figure
boxplot(RMSE_PDS,name)
ylabel('RMSE')
set(gca, 'FontSize', 16)



RMSE=[RMSE_DS,RMSE_PDS];
name=["ABIAN-DS","CARMONA-DS","SISAV-DS","ABIAN-PDS","CARMONA-PDS","SISAV-PDS"];
figure
boxplot(RMSE,name)
ylabel('RMSE')
set(gca, 'FontSize', 16)



T_ABIAN_R_PDS   =   Yhat_ABIAN_PDS(:,index_CeC)     * model_PCA.P(:,1:model_PCA.k);
T_CARMONA_R_PDS =   Yhat_CARMONA_PDS(:,index_CeC)   * model_PCA.P(:,1:model_PCA.k);
% T_CeC_R_DS     =   x_CeC_transfer(:,index_CeC) * model_PCA.P(:,1:model_PCA.k);
T_SISAV_R_PDS   =   Yhat_SISAV_PDS(:,index_CeC)     * model_PCA.P(:,1:model_PCA.k);

figure
hold on
plot(T_ABIAN_R_PDS(:,1),T_ABIAN_R_PDS(:,2),'b.')
plot(T_CARMONA_R_PDS(:,1),T_CARMONA_R_PDS(:,2),'r*')
plot(T_CeC_1(:,1),T_CeC_1(:,2),'ks')
plot(T_SISAV_R_PDS(:,1),T_SISAV_R_PDS(:,2),'go')
hold off
xlabel('PC1')
ylabel('PC2')
legend('ABIAN','CARMONA','CeC','SISAV','Location','southeast')
axis([-2000 1000 -1000  500 ])
set(gca, 'FontSize', 16)
box on



figure
subplot(2,2,1)
hold on
title('ABIAN')
plot(T_CeC_1(:,1),T_CeC_1(:,2),'ks')
plot(T_ABIAN_1(:,1),T_ABIAN_1(:,2),'b.')
plot(T_ABIAN_R_DS(:,1),T_ABIAN_R_DS(:,2),'r*')
plot(T_ABIAN_R_PDS(:,1),T_ABIAN_R_PDS(:,2),'go')
hold off
xlabel('PC1')
ylabel('PC2')
axis([-2000 1000 -1000  500 ])
set(gca, 'FontSize', 16)
box on
subplot(2,2,2)
hold on
title('CARMONA')
plot(T_CeC_1(:,1),T_CeC_1(:,2),'ks')
plot(T_CARMONA_1(:,1),T_CARMONA_1(:,2),'b.')
plot(T_CARMONA_R_DS(:,1),T_CARMONA_R_DS(:,2),'r*')
plot(T_CARMONA_R_PDS(:,1),T_CARMONA_R_PDS(:,2),'go')
hold off
xlabel('PC1')
ylabel('PC2')
axis([-2000 1000 -1000  500 ])
set(gca, 'FontSize', 16)
box on
subplot(2,2,3)
hold on
title('SISAV')
plot(T_CeC_1(:,1),T_CeC_1(:,2),'ks')
plot(T_SISAV_1(:,1),T_SISAV_1(:,2),'b.')
plot(T_SISAV_R_DS(:,1),T_SISAV_R_DS(:,2),'r*')
plot(T_SISAV_R_PDS(:,1),T_SISAV_R_PDS(:,2),'go')
hold off
xlabel('PC1')
ylabel('PC2')
axis([-2000 1000 -1000  500 ])
set(gca, 'FontSize', 16)
box on

subplot(2,2,4)
hold on
plot(T_CeC_1(1,1),T_CeC_1(1,2),'ks')
plot(T_SISAV_1(1,1),T_SISAV_1(1,2),'b.')
plot(T_SISAV_R_DS(1,1),T_SISAV_R_DS(1,2),'r*')
plot(T_SISAV_R_PDS(1,1),T_SISAV_R_PDS(1,2),'go')
legend('reference','initial','DS','PDS','Location','southeast')
set(gca, 'FontSize', 16)
box on
