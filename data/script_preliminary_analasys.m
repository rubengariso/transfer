clc; close all; clear; tic

%% Preambulo

addpath ../matalab_KEYS_funcion/'PCA Package Code'
addpath ../matalab_KEYS_funcion/functions/

%% Ler & organizar dados

raw_data=load('data/Coagulacoes_FTIR_C&C_24/C&C_2024_11_11.mat','CeC');
CeC_raw.lambda=raw_data.CeC.lambda*10^(-2);
CeC_raw.sample=raw_data.CeC.name;
CeC_raw.spectra=10^-2*raw_data.CeC.x;

raw_data=load('data/ABIAM/Resultados FTIR/ABIAN.mat','ABIAN');
ABIAN_raw.lambda=raw_data.ABIAN.lambda;
ABIAN_raw.sample=raw_data.ABIAN.sample;
ABIAN_raw.spectra=raw_data.ABIAN.spectra;

raw_data=load('data/Carmona/carmona.mat','CARMONA');
CARMONA_raw.lambda=raw_data.CARMONA.lambda;
CARMONA_raw.sample=raw_data.CARMONA.sample;
CARMONA_raw.spectra=raw_data.CARMONA.spectrum;

load('data/SISAV/SISAV.mat','x','NAME','lambda')
SISAV_raw.spectra=100*x;
SISAV_raw.sample=categorical(NAME);
SISAV_raw.lambda=lambda*10^(-6);

clear x NAME lambda raw_data

numberOfReplicates=3;

%% Uniformizar label das amostras

load('data/chave_das_amostras.mat')
sample_label_match=tranfer_v2;


for i=1:size(sample_label_match,1)

    idx_A=ismember(ABIAN_raw.sample, sample_label_match(i,1));
    ABIAN_raw.sample(idx_A)=sample_label_match(i,2);

    idx_A=ismember(CARMONA_raw.sample, sample_label_match(i,1));
    CARMONA_raw.sample(idx_A)=sample_label_match(i,2);

    idx_A=ismember(CeC_raw.sample, sample_label_match(i,1));
    CeC_raw.sample(idx_A)=sample_label_match(i,2);

    idx_A=ismember(SISAV_raw.sample, sample_label_match(i,1));
    SISAV_raw.sample(idx_A)=sample_label_match(i,2);

end

%% reorganizar os espectro de modo a ter as amostras comuns

ID_global=[SISAV_raw.sample; ABIAN_raw.sample; CARMONA_raw.sample; CeC_raw.sample];
ID_unique=unique(ID_global);

idx_A=ismember(ID_unique,ABIAN_raw.sample); idx_B=ismember(ID_unique,CARMONA_raw.sample);
idx_C=ismember(ID_unique,CeC_raw.sample);     idx_D=ismember(ID_unique,SISAV_raw.sample);

idx=[idx_A, idx_B, idx_C, idx_D];

ID_comum=ID_unique(sum(idx,2)==4);

numberOfComumSamples=size(ID_comum,1);


lambda.ABIAN=ABIAN_raw.lambda(:);
lambda.CARMONA=CARMONA_raw.lambda(:);
lambda.CeC=CeC_raw.lambda(:);
lambda.SISAV=SISAV_raw.lambda(:);

spectra.ABIAN=zeros(numberOfReplicates*numberOfComumSamples,size(ABIAN_raw.spectra,2));
spectra.CARMONA=zeros(numberOfReplicates*numberOfComumSamples,size(CARMONA_raw.spectra,2));
spectra.CeC=zeros(numberOfReplicates*numberOfComumSamples,size(CeC_raw.spectra,2));
spectra.SISAV=zeros(numberOfReplicates*numberOfComumSamples,size(SISAV_raw.spectra,2));



for i=1:numberOfComumSamples
    index=(i-1)*numberOfReplicates+1:i*numberOfReplicates;

    % ABIAN
    idx_A=ismember(ABIAN_raw.sample, ID_comum(i));
    spectra.ABIAN(index,:)=ABIAN_raw.spectra(idx_A,:);

    % CARMONA
    idx_A=ismember(CARMONA_raw.sample, ID_comum(i));
    spectra.CARMONA(index,:)=CARMONA_raw.spectra(idx_A,:);
    
    % CORREIA&CORREIA
    idx_A=ismember(CeC_raw.sample, ID_comum(i));
    spectra.CeC(index,:)=CeC_raw.spectra(idx_A,:);

    % SISAV
    idx_A=ismember(SISAV_raw.sample, ID_comum(i));
    spectra.SISAV(index,:)=SISAV_raw.spectra(idx_A,:);

end

%% Modelo PCA

% definir a gama de lambdas a usar

lambda_inicio=min([CARMONA_raw.lambda(1);SISAV_raw.lambda(1);ABIAN_raw.lambda(1);CeC_raw.lambda(1)]);
lambda_fim=max([CARMONA_raw.lambda(end);SISAV_raw.lambda(end);ABIAN_raw.lambda(end);CeC_raw.lambda(end)]);

index_CeC=find(CeC_raw.lambda<=lambda_inicio & CeC_raw.lambda>=lambda_fim);
lambda_ref=CeC_raw.lambda(index_CeC)';

numberOfVariables=size(lambda_ref,2);

% espectros estimados usando splines

spectra_splines.CARMONA   = spline(lambda.CARMONA,spectra.CARMONA,lambda_ref);
spectra_splines.SISAV     = spline(lambda.SISAV,spectra.SISAV,lambda_ref);
spectra_splines.CeC       = spline(lambda.CeC,spectra.CeC,lambda_ref);
spectra_splines.ABIAN     = spline(lambda.ABIAN,spectra.ABIAN,lambda_ref);

figure
subplot(2,2,1)
plot(lambda_ref,spectra_splines.ABIAN)
xlabel('\lambda (cm^{-1})')
ylabel('T(%)')
title('ABIAN')
set(gca, 'FontSize', 16)
axis([500 4000 0 90])
subplot(2,2,2)
plot(lambda_ref,spectra_splines.CARMONA)
xlabel('\lambda (cm^{-1})')
ylabel('T(%)')
title('CARMONA')
axis([500 4000 0 90])
set(gca, 'FontSize', 16)
subplot(2,2,3)
plot(lambda_ref,spectra_splines.CeC)
xlabel('\lambda (cm^{-1})')
ylabel('T(%)')
title('C&C')
axis([500 4000 0 90])
set(gca, 'FontSize', 16)
subplot(2,2,4)
plot(lambda_ref,spectra_splines.SISAV)
xlabel('\lambda (cm^{-1})')
ylabel('T(%)')
title('SISAV')
set(gca, 'FontSize', 16)
axis([500 4000 0 90])

z_ABIAN=spectra_splines.ABIAN -mean(spectra_splines.ABIAN);
z_CARMONA=spectra_splines.CARMONA- mean(spectra_splines.CARMONA);
z_CeC=spectra_splines.CeC-mean(spectra_splines.CeC);
z_SISAV=spectra_splines.SISAV-mean(spectra_splines.SISAV);

x_pca=[z_ABIAN;z_CARMONA;z_CeC;z_SISAV];
alpha=1/100;
kSelect='kSelectCpv';
[ model_PCA ] = pcaModel(x_pca, kSelect, 'datamat' );


% Determinar PC para cada laboratorio
T_ABIAN_1   =   z_ABIAN     * model_PCA.P(:,1:model_PCA.k);
T_CARMONA_1 =   z_CARMONA   * model_PCA.P(:,1:model_PCA.k);
T_CeC_1     =   z_CeC       * model_PCA.P(:,1:model_PCA.k);
T_SISAV_1   =   z_SISAV     * model_PCA.P(:,1:model_PCA.k);


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
% axis([-2000 1000 -1000  500 ])
box on

%% Guardar outputs
% save('data/dados_organizados.mat','spectra','lambda','ID_comum');
% save('data/pca_model.mat','model_PCA','index_CeC','T_ABIAN_1','T_CARMONA_1','T_CeC_1','T_SISAV_1');

