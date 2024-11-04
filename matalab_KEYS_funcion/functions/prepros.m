function [ppXtrain,ppXtest] = prepros(Xtrain,Xtest,pp_metodo,sgd_model)
%UNTITLED2 Summary of this function goes here

%INPUTS
%  Xtest - conjunto de teste
%  Xtest - conjunto de teste
%  pp_metodo tecnica - qual o metodo que se pertende
% parametros '
%   Detailed explanation goes here´


switch pp_metodo
    case 'mean_centering'
        mc_model.mean=mean(Xtrain);
        ppXtrain=Xtrain-mc_model.mean;
        %test
        ppXtest=Xtest-mc_model.mean;
    case 'MSC'
        msc_model.first=1; % first: first variable used for correction
        msc_model.last=size(Xtrain,2);% last: last variable used for correction
        [Xpp,msc_model.me]=msc(Xtrain,msc_model.first,msc_model.last);
        msc_model.mean=mean(Xpp);
        ppXtrain=Xpp-msc_model.mean;

    % test
        [~,~,Xpp]=msc(msc_model.me,msc_model.first,msc_model.last, Xtest);
        ppXtest=Xpp-msc_model.mean;       
    case 'SNV'
        [Xpp]=snv(Xtrain);
        snv_model.mean=mean(Xpp);
        ppXtrain=Xpp-snv_model.mean;
        [Xpp]=snv(Xtest);
        ppXtest=Xpp-snv_model.mean; 
    case 'SGV'
%         [Xp]=snv(Xtrain);
%         snv_model.mean=mean(Xp);
        ppXtrain_0=Xtrain; % or MSC or SNV
        Xpp= deriv(ppXtrain_0,sgd_model.der,sgd_model.window,sgd_model.order);
        sgd_model.mean=mean(Xpp);
        ppXtrain=Xpp-sgd_model.mean; 
        
        % test
        ppXtest_0=Xtest; % or MSC or SNV (corresponding to train)
        Xpp = deriv(ppXtest_0,sgd_model.der,sgd_model.window,sgd_model.order);
        ppXtest=Xpp-sgd_model.mean;
        case 'SNV-SGV'
        [Xpp]=snv(Xtrain);
%         snv_model.mean=mean(Xpp);
        ppXtrain_0=Xpp;
        [Xpp]=snv(Xtest);
        ppXtest_0=Xpp; 
%         ppXtrain_0=Xtrain; % or MSC or SNV
        Xpp= deriv(ppXtrain_0,sgd_model.der,sgd_model.window,sgd_model.order);
        sgd_model.mean=mean(Xpp);
        ppXtrain=Xpp-sgd_model.mean; 
        
        Xpp = deriv(ppXtest_0,sgd_model.der,sgd_model.window,sgd_model.order);
        ppXtest=Xpp-sgd_model.mean;
    otherwise
        fprintf('Wrong data pretreat method!')
end






end

