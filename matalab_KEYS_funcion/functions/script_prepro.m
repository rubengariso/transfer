% Script to exemplify the use of multiple pre-processing methos


clear
clc
close all
 
%% Load Data ==============================================================

load spectra.mat
Xtrain=NIR(1:40,:);
Ytrain=octane(1:40);

Xtest=NIR(41:end,:);
Ytest=octane(41:end);

%==========================================================================

%% Mean-centering =========================================================

%train
mc_model.mean=mean(Xtrain);
ppXtrain=Xtrain-mc_model.mean;

%test
ppXtest=Xtest-mc_model.mean;
            
%==========================================================================

%% Multiplicative scatter correction: MSC =================================

%train
% msc followed by mean-centering
msc_model.first=1; % first: first variable used for correction
msc_model.last=size(Xtrain,2);% last: last variable used for correction
[Xpp,msc_model.me]=msc(Xtrain,msc_model.first,msc_model.last);
msc_model.mean=mean(Xpp);
ppXtrain=Xpp-msc_model.mean;

% test
[~,~,Xpp]=msc(msc_model.me,msc_model.first,msc_model.last, Xtest);
ppXtest=Xpp-msc_model.mean;
            
%==========================================================================

%% Standard normal variate: SNV ===========================================

%train
% snv followed by mean-centering
[Xpp]=snv(Xtrain);
snv_model.mean=mean(Xpp);
ppXtrain=Xpp-snv_model.mean; 

% test
[Xpp]=snv(Xtest);
ppXtest=Xpp-snv_model.mean; 
        
%==========================================================================

%% Savitzky-Golay differentiation: SG =====================================

% train
% SGV followed by mean-centering
ppXtrain_0=Xtrain; % or MSC or SNV
sgd_model.der=1 ;%degree of the derivative; it must be <= order
sgd_model.window=7;% (optional), (1x1) the number of points in filter, it must be >3 and odd
sgd_model.order=1;%  (optional), (1x1) the order of the polynomial. It must be <=5 and <= (window-1)
Xpp= deriv(ppXtrain_0,sgd_model.der,sgd_model.window,sgd_model.order);
sgd_model.mean=mean(Xpp);
ppXtrain=Xpp-snv_model.mean; 

% test
ppXtest_0=Xtest; % or MSC or SNV (corresponding to train)
Xpp = deriv(ppXtest_0,sgd_model.der,sgd_model.window,sgd_model.order);
ppXtest=Xpp-sgd_model.mean; 

%==========================================================================