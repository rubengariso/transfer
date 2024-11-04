function [model] = rpcaModel(Xt, Xv, aG, varargin)
%RPCAMODEL outputs an RPCA model for subsequent online monitoring. This
%function can also perform parameter selection for the the alpha and eta
%parameters. The user can also specify the parameter values instead of
%automatically selecting them.
%
% References: 
%     Li, W., Yue, H. H., Valle-Cervantes, S., Qin, S. J.,
%     "Recursive PCA for adaptive process monitoring",
%     Journal of Process Control, 5, 471-486

%
% Required input arguments: 
%   Xt   : A data matrix of observations in which rows represent observations, and columns represent variables.
%          Xt is used to obtain an initial PCA model using pcaModel.
%   Xv   : A data matrix used to validate parameter estimates.
%   aG   : The desired global false detection rate.
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
% I/O: model=rpcaModel( Xt, Xv, aG, 'eta', [], 'aT', [],'aE', []);
%
% The output of RPCAMODEL is a structure containing:
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
%   model.uclT        : The upper control limit of the T^2-statistic.
%   model.uclE        : The upper control limit of the Q-statistic.
%   model.aG          : The desired global false detection rate
%   model.T           : The current score vector.
%   model.loc_old     : The sample mean from the previous in-control period.
%
% see also pcaModel, rpcaMonitor, aSolve, fSpeedSolve.

aT = NaN;
aE = NaN;
model.eta = []; 
cpvThresh = 0.95;
typeLimit = 'Theoretic';
if nargin>3,
    varargin=reshape(varargin,2,[]);
    for c=1:size(varargin,2),
        switch varargin{1,c}
            case 'eta'
                model.eta = varargin{2,c};
            case 'aT'
                aT=varargin{2,c};
            case 'aE'
                aE=varargin{2,c};
            case 'cpvThresh'
                cpvThresh=varargin{2,c};
            case 'typeLimit'
                typeLimit=varargin{2,c};
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PART I
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    
    %model.eta = eta;

%     n = size(Xt,1);
%     model.n = n;
%     model.N = n;
%     %standardize the data
%     model.loc = mean(Xt);
%     model.R = corr(Xt);
%     model.S = cov(Xt);
    model.kSelectType = 'covmat';
    model.kSelect = 'kSelectCpv';
    model.type = 'RPCA';
    model.cpvThresh = cpvThresh;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PART II
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Fit the initial PCA model and set control limits.
%fSpeedSolve will compute an initial RPCA model for each forgetting factor
    if isempty(model.eta) == 1, 
    [model.eta] = fSpeedSolve(Xt,Xv, model);
    end
%Here, we fit the initial weighted PCA model with the corrrect forgetting
%factor.
    model = wPcaModel(Xt, model.eta, model.kSelect, model.kSelectType, 'cpvThresh', cpvThresh);
    model.aT = aT;
    model.aE = aE;
    model.aG = aG;
    %%%%%
%     model.aT = eps;
%     model.aE = eps;
%     model.ucl.sT = Inf;
%     model.ucl.sE = Inf;
%     [STAT, UCL, model] = rpcaMonitor(Xv, model);

    switch typeLimit
        case 'Theoretic'
            model.aT = aG/2;
            model.aE = aG/2;
            model.ucl.sT = limitT(model.aT,model.k);
            model.ucl.sE = limitQ(model.L, model.k, model.aE);
        case 'Empiric'
            %Find optimal alpha values
            [model] = aSolve(Xv, model, aG, model.aT, model.aE);
    end
    %Estimate final model using optimal values of alpha
    [~, ~, model] = rpcaMonitor(Xv, model);
end

