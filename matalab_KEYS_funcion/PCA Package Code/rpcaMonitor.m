function [stat, ucl, model] = rpcaMonitor(X, model)
%RPCAMONITOR monitors process data and updates the RPCA model to the
%evolution of the process. Updated model parameters, monitoring statistics,
%and control limits are returned as output.
%
% References: 
%     Li, W., Yue, H. H., Valle-Cervantes, S., Qin, S. J.,
%     "Recursive PCA for adaptive process monitoring",
%     Journal of Process Control, 5, 471-486
%
% Required input arguments: 
%   X     : A data matrix of observations in which rows represent observations, and columns represent variables.
%          X are the observations to be monitored.
%   model : An 'RPCA' model structure to start monitoring from.
%
% Optional input arguments:
%   eta : The user provided value of the forgetting factor eta.
%            In principle, eta can values such that 0 <= eta <= 1.
%            Though a practical lower bound is 0.99.
%   aT     : The user provided value of the T^2-statistic alpha
%            parameter. aT takes values such that 0 <= aT <= 1.
%   aE     : The user provided value of the model eQ-statistic alpha
%            parameter. aE takes values such that 0 <= aT <= 1.
%
% I/O: [stat, ucl, model]=rpcaModel(X, model);
%
% The output of RPCAMONITOR is a structure containing:
%   model.type        : The type of model being computed is stored as 'RPCA'.
%   model.aT          : Value of the T^2-statistic alpha parameter.
%   model.aE          : Value of the Q-statistic alpha parameter.
%   model.eta      : Value of eta, the forgetting parameter.
%   model.n           : The number of observations evaluated by the RPCA
%                       model up to the current time.
%   model.N           : The size of the Xt initially evaluated by pcaModel.
%   model.loc         : The weighted sample mean.
%   model.R           : The weighted sample correlation matrix.
%   model.S           : The weighted sample covariance matrix.
%   model.kSelectType : Latent variable selection always performed using
%                       'covmat' as the input type.
%   model.kSelect     : Latent variable selection always performed using
%                       'CPV' as the seelction method.
%   model.P           : The PCA model loadings.
%   model.L           : The singular values.
%   model.k           : Latent variable selection always performed using
%                       'covmat'.
%   model.ucl.sT        : The upper control limit of the T^2-statistic.
%   model.ucl.sE        : The upper control limit of the Q-statistic.
%   model.aG          : The desired global false detection rate
%   model.T           : The current score vector.
%   model.loc_old     : The sample mean from the previous in-control period.
%   model.kVec        : Vector of old number of components retained.
%   stat.sT           : Vector of the T^2 statistics of the observations in X.
%   stat.sE           : Vector of the Q-statistics of the observations in X.
%   ucl.ucl.sT          : Vector of the T^2-statistic control limits used to monitor the observations in X. 
%   ucl.ucl.sE          : Vector of the Q-statistic control limits used to monitor the observations in X. 
%
% see also pcaModel, rpcaModel, aSolve, fSpeedSolve.

    n = size(X,1);
    stat.sE = NaN(n,1);
    stat.sT =  NaN(n,1);
    ucl.sE =  NaN(n,1);
    ucl.sT =  NaN(n,1);
    model.kVec = [];
    for i = 1:n,
        eta_n = model.n/(model.n+1)*model.eta; % weight parameter
        xi = X(i,:); % in LV, this will be the signal from the new seal
        xi_center = xi - model.loc;% mean centering of the new signal
        xi_scale = (xi_center)./sqrt(diag(model.S))'; 
        model.T =  xi_scale*model.P(:,1:model.k); % PCA score of the new signal
        sE = sum(( xi_scale -  xi_scale*model.P(:,1:model.k)*model.P(:,1:model.k)')'.^2); % Q statistic, in LV, the index i for res is not needed, we will show res in a waveform chart
        sT = sum(((xi_scale* model.P(:,1: model.k)*inv(model.L(1: model.k,1: model.k).^.5)).^2),2);
        stat.sE(i) = sE;
        stat.sT(i) = sT;

        if sE<=model.ucl.sE & sT<=model.ucl.sT, 
            model.loc_old = model.loc; % store the old mean vector b
            model.loc = eta_n*model.loc + (1-eta_n)*xi; % update b (this should be done using the raw data, not a centered version)
            d_loc = model.loc - model.loc_old; % difference between new and old b
            model.S = eta_n*(model.S+d_loc'*d_loc)+(1-eta_n)*xi_center'*xi_center; % update for the covariance matrix to be applied for scaling
            model.R = corrcov(model.S);
            model.n = model.n + 1; % update of the number of samples
            model_i = pcaModel(model.R, model.kSelect, model.kSelectType, 'cpvThresh', model.cpvThresh);
            model.P = model_i.P;
            model.L = model_i.L;
            model.k = model_i.k;
           
% model.k = 5; %changed for test
            model.ucl.sT = limitT(model.aT,model.k);
            model.ucl.sE = limitQ(model.L, model.k, model.aE);
        end
        ucl.sE(i) = model.ucl.sE;
        ucl.sT(i) = model.ucl.sT;
        model.kVec = [model.kVec, model.k];
    end
end