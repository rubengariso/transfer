clc;clear; close all;tic
load("transfer_calibration.mat")

addpath ../matalab_KEYS_funcion/functions/
addpath ../matalab_KEYS_funcion/'PCA Package Code'/

%% Data augmentation

x_CeC_augmentation=zeros(3*size(x_CeC_transfer,1),size(x_CeC_transfer,2));
x_ABIAN_augmentation=zeros(3*size(x_ABIAN_transfer,1),size(x_ABIAN_transfer,2));
x_CARMONA_augmentation=zeros(3*size(x_CARMONA_transfer,1),size(x_CARMONA_transfer,2));
x_SISAV_augmentation=zeros(3*size(x_SISAV_transfer,1),size(x_SISAV_transfer,2));
IDX_start=1:3:3*size(x_SISAV_transfer,1);
iiii=1:size(x_SISAV_transfer,1);
IDX_start_2=1:9:3*size(x_CeC_transfer,1);

for i=1:size(x_SISAV_transfer,1)
    first_Cec = x_CeC_transfer(i, :);
    x_CeC_augmentation(IDX_start(i):IDX_start(i)+2,:)   =   repmat(first_Cec, 3, 1);
    iiii_CeC(IDX_start(i):IDX_start(i)+2,:)   =   repmat(iiii(i), 3, 1);
end

for i=1:size(IDX_start_2,2)
    first_ABIAN     = x_ABIAN_transfer(IDX_start(i):IDX_start(i)+2, :);
    first_CARMONA   = x_CARMONA_transfer(IDX_start(i):IDX_start(i)+2, :);
    first_SISAV     = x_SISAV_transfer(IDX_start(i):IDX_start(i)+2, :);
    frist_iii       = iiii(IDX_start(i):IDX_start(i)+2);
    id_replica      = ID_final_dos_finais(IDX_start(i):IDX_start(i)+2, :);
    x_ABIAN_augmentation(IDX_start_2(i):IDX_start_2(i)+8,:)     =  repmat(first_ABIAN, 3, 1);
    x_CARMONA_augmentation(IDX_start_2(i):IDX_start_2(i)+8,:)   =  repmat(first_CARMONA, 3, 1);
    x_SISAV_augmentation(IDX_start_2(i):IDX_start_2(i)+8,:)     =  repmat(first_SISAV, 3, 1);
    iiii_ABIAN(IDX_start_2(i):IDX_start_2(i)+8,:)               =   repmat(frist_iii', 3, 1);
    ID(IDX_start_2(i):IDX_start_2(i)+8,:)                       =   repmat(id_replica, 3, 1);
end


n_random=200;
ID_unique=unique(ID);
matrix=zeros(size(ID_unique,1),n_random);

for col = 1:n_random
    random_indices = randperm(size(ID_unique,1),round(0.2*size(ID_unique,1)));
    matrix(random_indices, col) = 1;  % Set those indices to 1 in the column
end
matrix_end=nan(size(ID,1),n_random);
for i=1:size(matrix,1)
    vector=matrix(i,:);
     matrix_end(IDX_start_2(i):IDX_start_2(i)+8,:)=repmat(vector, 9, 1);
end
toc

sum(sum(isnan(matrix_end))) 
% if matrix==0 => calibration
% if matrix==1 => validation

[row, col] = find(matrix_end == 1);
n = ceil(length(row) /length(unique(col)));
validation = NaN(n, 200);

for i = 1:200
    col_indices = find(col == i);
    validation(1:n, i) = row(col_indices);
end

[row, col] = find(matrix_end == 0);
n = ceil(length(row) /length(unique(col)));
calibration = NaN(n, 200);

for i = 1:200
    col_indices = find(col == i);
    calibration(1:n, i) = row(col_indices);
end

index_mccv.train=calibration;
index_mccv.test=validation;
%% DS
toc
[F_DS_ABIAN,model_ABIAN,RMSE_ABIAN,Yhat_ABIAN]          = direct_standardization(x_CeC_augmentation,x_ABIAN_augmentation,index_mccv);
toc
[F_DS_CARMONA,model_CARMONA,RMSE_CARMONA,Yhat_CARMONA]  = direct_standardization(x_CeC_augmentation,x_CARMONA_augmentation,index_mccv);
toc
[F_DS_SISAV,model_SISAV,RMSE_SISAV,Yhat_SISAV]          = direct_standardization(x_CeC_augmentation,x_SISAV_augmentation,index_mccv);
toc

RMSE_DS=[RMSE_ABIAN, RMSE_CARMONA, RMSE_SISAV];

name=["ABIAN-DS","CARMONA-DS","SISAV-DS"];
figure
boxplot(RMSE_DS,name)
ylabel('RMSE')
set(gca, 'FontSize', 16)



T_ABIAN_R   =   Yhat_ABIAN(:,index_CeC)     * model_PCA.P(:,1:model_PCA.k);
T_CARMONA_R =   Yhat_CARMONA(:,index_CeC)   * model_PCA.P(:,1:model_PCA.k);
T_CeC_R     =   x_CeC_transfer(:,index_CeC) * model_PCA.P(:,1:model_PCA.k);
T_SISAV_R   =   Yhat_SISAV(:,index_CeC)     * model_PCA.P(:,1:model_PCA.k);

figure
hold on
plot(T_ABIAN_1(:,1),T_ABIAN_1(:,2),'b.')
plot(T_CARMONA_R(:,1),T_CARMONA_R(:,2),'r*')
plot(T_CeC_R(:,1),T_CeC_R(:,2),'ks')
plot(T_SISAV_R(:,1),T_SISAV_R(:,2),'go')
hold off
xlabel('PC1')
ylabel('PC2')
legend('ABIAN','CARMONA','CeC','SISAV')
set(gca, 'FontSize', 16)
box on
% 

%% PCS
window=10;
toc
[F_ABIAN_PDS,RMSE_ABIAN_PDS,Yhat_ABIAN_PDS] = piecewise_direct_standardization_v2(x_ABIAN_augmentation,x_CeC_augmentation,window,CeC.lambda,ABIAN.lambda,index_mccv);
toc
[F_CARMONA_PDS,RMSE_CARMONA_PDS,Yhat_CARMONA_PDS] = piecewise_direct_standardization_v2(x_CARMONA_augmentation,x_CeC_augmentation,window,CeC.lambda,CARMONA.lambda,index_mccv);
toc
[F_SISAV_PDS,RMSE_SISAV_PDS,Yhat_SISAV_PDS] = piecewise_direct_standardization_v2(x_SISAV_augmentation,x_CeC_augmentation,window,CeC.lambda,SISAV.lambda,index_mccv);
toc
RMSE_PSD=[RMSE_SISAV_PDS,RMSE_CARMONA_PDS,RMSE_ABIAN_PDS];
name=["SISAV-PDS","CARMONA-PDS","ABIAN-PDS"];
figure
boxplot(RMSE,name)
ylabel('RMSE')
set(gca, 'FontSize', 16)



RME=[ RMSE_DS,RMSE_PSD];
name=["ABIAN-DS","CARMONA-DS","SISAV-DS","SISAV-PDS","CARMONA-PDS","ABIAN-PDS"];
figure
boxplot(RME,name)
ylabel('RMSE')
set(gca, 'FontSize', 16)
toc
