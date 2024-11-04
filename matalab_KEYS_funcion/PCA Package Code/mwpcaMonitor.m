function [stat, ucl, model] = mwpcaMonitor(X, model)
%MWPCAMODEL monitors process data and updates the MWPCA model to the
%evolution of the process. Updated model parameters, monitoring statistics,
%and control limits are returned as output.
%
% References: 
%     Wang, X., Kruger, U., Irwin, G. W.,
%     "Process Monitoring Approach Using Fast Moving Window PCA",
%     Industrial & Engineering Chemistry Research, 15, 5691-5702
%
% Required input arguments: 
%   X     : A data matrix of observations in which rows represent observations, and columns represent variables.
%          X are the observations to be monitored.
%   model : An 'MWPCA' model structure to start monitoring from.
%
% I/O: [stat, ucl, model]=mwpcaMonitor(X, model);
%
% The output of MWPCAMONITOR is a structure containing:
%   model.type        : The type of model being computed is stored as 'MWPCA'.
%   model.aT          : Value of the T^2-statistic alpha parameter.
%   model.aE          : Value of the Q-statistic alpha parameter.
%   model.h           : Value of h, the window size parameter.
%   model.n           : The number of observations from the dataset Xt used
%                       to train the initial PCA model. Typically, this equals h, but may be
%                       less if a selected h is larger than the number of observations in Xt.
%   model.loc         : The weighted sample mean.
%   model.R           : The weighted sample correlation matrix.
%   model.sd          : The weighted sample standard deviations.
%   model.kSelectType : Latent variable selection always performed using
%                       'covmat' as the input type.
%   model.kSelect     : Latent variable selection always performed using
%                       'CPV' as the seelction method.
%   model.P           : The PCA model loadings.
%   model.L           : The singular values.
%   model.k           : Latent variable selection always performed using 'covmat'
%   model.kVec        : Vector of old number of components retained.
%   model.ucl.sT      : The upper control limit of the T^2-statistic.
%   model.ucl.sE      : The upper control limit of the Q-statistic.
%   model.aG          : The desired global false detection rate
%   model.T           : The current score vector.
%   model.Xh          : The observations in the current, updated window.
%   model.Xh_old      : Used to store observations from the previous period
%                       window.
%   stat.sT           : Vector of the T^2 statistics of the observations in X.
%   stat.sE           : Vector of the Q-statistics of the observations in X.
%   ucl.sT            : Vector of the T^2-statistic control limits used to monitor the observations in X. 
%   ucl.sE            : Vector of the Q-statistic control limits used to monitor the observations in X. 
%
% see also pcaModel, mwpcaModel, aSolve, fSpeedSolve.


    n = size(X,1);
    stat.sE = NaN(n,1);
    stat.sT =  NaN(n,1);
    ucl.sE =  NaN(n,1);
    ucl.sT =  NaN(n,1);
    model.kVec = [];
    for i = 1:n,
        xi = X(i,:); % the new observation
        xi_scale = (xi - model.loc)./model.sd; 
        model.T =  xi_scale*model.P(:,1:model.k); % PCA score of the new signal
        sE = sum(( xi_scale -  xi_scale*model.P(:,1:model.k)*model.P(:,1:model.k)')'.^2); % Q statistic, in LV, the index i for res is not needed, we will show res in a waveform chart
        sT = sum(((xi_scale* model.P(:,1: model.k)*inv(model.L(1: model.k,1: model.k).^.5)).^2),2);
        stat.sE(i) = sE;
        stat.sT(i) = sT;     
        
        if sE<=(model.ucl.sE) & sT<=model.ucl.sT, 
            model.Xh = [model.Xh_old(2:end,:);xi];
            
            [loc, sd, R ] = mwcov( xi, model.Xh_old(1,:), model.loc, model.sd, model.R, model.h);
            model.loc = loc';
            model.sd = sd'; 
            model.R = R';
            
%             model.sd = std(model.Xh);
%             model.R = corr(model.Xh);
%             model.loc = mean(model.Xh);
                        
            if model.n < model.h,
            model.n = model.n + 1; % update of the number of samples
            end
            model_i = pcaModel(model.R, model.kSelect, model.kSelectType, 'cpvThresh', model.cpvThresh);
            model.P = model_i.P;
            model.L = model_i.L;
            model.k = model_i.k;
            model.ucl.sT = limitT(model.aT,model.k);
            model.ucl.sE = limitQ(model.L, model.k, model.aE);
            model.Xh_old = [model.Xh_old(2:end,:);xi];
        end
        ucl.sE(i) = model.ucl.sE;
        ucl.sT(i) = model.ucl.sT;
        model.kVec = [model.kVec, model.k];
    end
end


function [ loc_new, s_new, R_new ] = mwcov( x_new, x_old, loc_old, s_old, R_old, H)
% 24.10.2012

% Wang (2005)
    x_new=reshape(x_new,[],1);
    x_old=reshape(x_old,[],1);
    loc_old=reshape(loc_old,[],1);
    s_old=reshape(s_old,[],1);

    S_old=diag(s_old);
    iS_old=diag(1./s_old);
    b_til=(H*loc_old-x_old)/(H-1);%(1)
    Db_til=loc_old-b_til;%(2)
    x_k=iS_old*(x_old-loc_old);%(3)
    A=R_old-iS_old*(Db_til*Db_til')*iS_old-(x_k*x_k')/(H-1);%(4)
    loc_new=((H-1)*b_til+x_new)/H;%(5)
    Db=loc_new-b_til;%(6)
    s2=s_old.^2+Db.^2-Db_til.^2+((x_new-loc_new).^2-(x_old-loc_old).^2)/(H-1);%(7)
    s_new=sqrt(s2);
    % S_new=diag(s_new);%(8)
    iS_new=diag(1./s_new);
    x_L=iS_new*(x_new-loc_new);%(9)
    R_new=iS_new*S_old*A*S_old*iS_new+iS_new*(Db*Db')*iS_new+(x_L*x_L')/(H-1);%(10)


end

