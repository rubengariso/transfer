clc;clear; close all;tic

%% Preambulo

addpath ../matalab_KEYS_funcion/functions/
addpath ../matalab_KEYS_funcion/'PCA Package Code'/

% Ler dados

load('data/dados_organizados.mat','spectra','lambda','ID_comum')
load('data/pca_model.mat','model_PCA','index_CeC','T_ABIAN_1','T_CARMONA_1','T_CeC_1','T_SISAV_1')

numberOfSamples=length(ID_comum);
numberOfReplicates=3;

%% Cruzamento das amostras

% mix_samples='no';
mix_samples='yes';

switch mix_samples
    case 'no'
        index_destination=(1:numberOfSamples*numberOfReplicates)';
        index_source=(1:numberOfSamples*numberOfReplicates)';
    case 'yes'

        A=repmat(1:numberOfReplicates:numberOfReplicates*numberOfSamples,numberOfReplicates,1);
        index_destination=[];
        for i=1:numberOfReplicates,
            index_destination=[index_destination;A-1+i];
        end
        index_destination=index_destination(:);

        index_source=repmat(reshape(1:numberOfSamples*numberOfReplicates,numberOfReplicates,[]),numberOfReplicates,1);
        index_source=index_source(:);
end


rng("default")

n_random=200;
fraction_calibration=0.80;
nCal=round(fraction_calibration*numberOfSamples);
calibration=false(length(index_destination),n_random);
for i=1:n_random
    random_indices = randperm(numberOfSamples,nCal);
    calibration(:,i)=sum(ceil(index_destination/numberOfReplicates)==random_indices,2)==1;
end

index_mccv.train=calibration;
index_mccv.test=~calibration;

%% DS ---------------------------------------------------------------------

name_labs={'ABIAN','CARMONA','SISAV'};

RMSE_bp_DS=nan(n_random,length(name_labs));
toc
for i=1:length(name_labs)
    X_source=spectra.(name_labs{i})(index_source,:);
    mew_X_source=mean(X_source);
    Z_source=X_source-mew_X_source;

    X_destination=spectra.CeC(index_destination,:);
    mew_X_destination=mean(X_destination);
    Z_destination=X_destination-mew_X_destination;

    [F_DS.(name_labs{i}),model_DS.(name_labs{i}),RMSE_DS.(name_labs{i}),Yhat_DS.(name_labs{i})] = direct_standardization(Z_source,Z_destination,index_mccv);
    model_DS.(name_labs{i}).mew_source=mew_X_source;
    model_DS.(name_labs{i}).mew_destination=mew_X_destination;

    Yhat_DS.(name_labs{i})=(Z_source*F_DS.(name_labs{i}))+mew_X_destination;

    RMSE_bp_DS(:,i)=RMSE_DS.(name_labs{i});
    toc
end


figure
boxplot(RMSE_bp_DS,string(name_labs')+'-DS')
ylabel('RMSE')
set(gca, 'FontSize', 16)

T_ABIAN_DS   =   (Yhat_DS.ABIAN(:,index_CeC)-mean(spectra.CeC(:,index_CeC)))     * model_PCA.P(:,1:model_PCA.k);
T_CARMONA_DS =   (Yhat_DS.CARMONA(:,index_CeC)-mean(spectra.CeC(:,index_CeC)))   * model_PCA.P(:,1:model_PCA.k);
T_CeC_DS     =   (spectra.CeC(:,index_CeC)-mean(spectra.CeC(:,index_CeC))) * model_PCA.P(:,1:model_PCA.k);
T_SISAV_DS   =   (Yhat_DS.SISAV(:,index_CeC)-mean(spectra.CeC(:,index_CeC)))     * model_PCA.P(:,1:model_PCA.k);

figure
hold on
plot(T_ABIAN_DS(:,1),T_ABIAN_DS(:,2),'b.')
plot(T_CARMONA_DS(:,1),T_CARMONA_DS(:,2),'r*')
plot(T_CeC_1(:,1),T_CeC_1(:,2),'ks')
plot(T_SISAV_DS(:,1),T_SISAV_DS(:,2),'go')
hold off
xlabel('PC1')
ylabel('PC2')
legend('ABIAN','CARMONA','CeC','SISAV')
set(gca, 'FontSize', 16)
box on

figure
for i=1:length(name_labs)
    subplot(1,3,i)
    plot(spectra.CeC','g')
    hold on
    plot(Yhat_DS.(name_labs{i})','r')
    hold off
    title(name_labs{i})
    set(gca, 'FontSize', 16)
    box on
end

%% PDS --------------------------------------------------------------------

window=20;
RMSE_bp_PDS=nan(n_random,length(name_labs));
toc
for i=1:length(name_labs)
    X_source=spectra.(name_labs{i})(index_source,:);
    mew_X_source=mean(X_source);
    Z_source=X_source-mew_X_source;
    lambda_source=lambda.(name_labs{i});

    X_destination=spectra.CeC(index_destination,:);
    mew_X_destination=mean(X_destination);
    Z_destination=X_destination-mew_X_destination;
    lambda_destination=lambda.CeC;

    [F_PDS.(name_labs{i}),RMSE_PDS.(name_labs{i}),Yhat_PDS.(name_labs{i})] = piecewise_direct_standardization(Z_source,Z_destination,lambda_source,lambda_destination,window,index_mccv);

    model_PDS.(name_labs{i}).mew_source=mew_X_source;
    model_PDS.(name_labs{i}).mew_destination=mew_X_destination;

    Yhat_PDS.(name_labs{i})=(Z_source*F_PDS.(name_labs{i}))+mew_X_destination;

    RMSE_bp_PDS(:,i)=RMSE_PDS.(name_labs{i});

    toc
end


figure
boxplot(RMSE_bp_PDS,string(name_labs')+'-PDS')
ylabel('RMSE')
set(gca, 'FontSize', 16)

T_ABIAN_PDS   =   (Yhat_PDS.ABIAN(:,index_CeC)-mean(spectra.CeC(:,index_CeC)))     * model_PCA.P(:,1:model_PCA.k);
T_CARMONA_PDS =   (Yhat_PDS.CARMONA(:,index_CeC)-mean(spectra.CeC(:,index_CeC)))   * model_PCA.P(:,1:model_PCA.k);
T_CeC_PDS     =   (spectra.CeC(:,index_CeC)-mean(spectra.CeC(:,index_CeC))) * model_PCA.P(:,1:model_PCA.k);
T_SISAV_PDS   =   (Yhat_PDS.SISAV(:,index_CeC)-mean(spectra.CeC(:,index_CeC)))     * model_PCA.P(:,1:model_PCA.k);

figure
hold on
plot(T_ABIAN_PDS(:,1),T_ABIAN_PDS(:,2),'b.')
plot(T_CARMONA_PDS(:,1),T_CARMONA_PDS(:,2),'r*')
plot(T_CeC_1(:,1),T_CeC_1(:,2),'ks')
plot(T_SISAV_PDS(:,1),T_SISAV_PDS(:,2),'go')
hold off
xlabel('PC1')
ylabel('PC2')
legend('ABIAN','CARMONA','CeC','SISAV')
set(gca, 'FontSize', 16)
box on

figure
for i=1:length(name_labs)
    subplot(1,3,i)
    plot(spectra.CeC','g')
    hold on
    plot(Yhat_PDS.(name_labs{i})','r')
    hold off
    title(name_labs{i})
    set(gca, 'FontSize', 16)
    box on
end

figure
boxplot([RMSE_bp_DS RMSE_bp_PDS],[string(name_labs')+'-DS';string(name_labs')+'-PDS'])
ylabel('RMSE')
set(gca, 'FontSize', 16)
toc

%% TSR --------------------------------------------------------------------

RMSE_bp_TSR=nan(n_random,length(name_labs));
toc
for i=1:length(name_labs)
    X_source=spectra.(name_labs{i})(index_source,:);
    mew_X_source=mean(X_source);
    Z_source=X_source-mew_X_source;
    lambda_source=lambda.(name_labs{i});

    X_destination=spectra.CeC(index_destination,:);
    mew_X_destination=mean(X_destination);
    Z_destination=X_destination-mew_X_destination;
    lambda_destination=lambda.CeC;

    [Xrec,RMSE_bp_TSR(:,i)]=pca_spectra_transfer(Z_source,Z_destination,index_mccv);
    Yhat_TSR.(name_labs{i})=Xrec+mew_X_destination;
    toc
end

T_ABIAN_TSR   =   (Yhat_TSR.ABIAN(:,index_CeC)-mean(spectra.CeC(:,index_CeC)))     * model_PCA.P(:,1:model_PCA.k);
T_CARMONA_TSR =   (Yhat_TSR.CARMONA(:,index_CeC)-mean(spectra.CeC(:,index_CeC)))   * model_PCA.P(:,1:model_PCA.k);
T_CeC_TSR     =   (spectra.CeC(:,index_CeC)-mean(spectra.CeC(:,index_CeC))) * model_PCA.P(:,1:model_PCA.k);
T_SISAV_TSR   =   (Yhat_TSR.SISAV(:,index_CeC)-mean(spectra.CeC(:,index_CeC)))     * model_PCA.P(:,1:model_PCA.k);

figure
hold on
plot(T_ABIAN_TSR(:,1),T_ABIAN_TSR(:,2),'b.')
plot(T_CARMONA_TSR(:,1),T_CARMONA_TSR(:,2),'r*')
plot(T_CeC_1(:,1),T_CeC_1(:,2),'ks')
plot(T_SISAV_TSR(:,1),T_SISAV_TSR(:,2),'go')
hold off
xlabel('PC1')
ylabel('PC2')
legend('ABIAN','CARMONA','CeC','SISAV')
set(gca, 'FontSize', 16)
box on


figure
for i=1:length(name_labs)
    subplot(1,3,i)
    plot(spectra.CeC','g')
    hold on
    plot(Yhat_TSR.(name_labs{i})','r')
    hold off
    title(name_labs{i})
    set(gca, 'FontSize', 16)
    box on
end

figure
rmse_labels=[string(name_labs')+'-DS';string(name_labs')+'-PDS';string(name_labs')+'-TSR'];
boxplot([RMSE_bp_DS RMSE_bp_PDS RMSE_bp_TSR],rmse_labels)
ylabel('RMSE')
set(gca, 'FontSize', 16)
toc
