function [model]=rplsModel(Xt, Xv, Yt,Yv,aG, k, varargin) % 
%PLSMODEL outputs an pls model for subsequent online monitoring. 
%
% References: 
%Add the algorithm paper...
%     Kourti, T., and MacGregor, J.,
%     "Process analysis, monitoring and diagnosis, using multivariate
%     projection methods",
%     Chemometrics and Intelligent Laboratory Systems, 28:1, pp. 3-21


%
% Required input arguments: 
%   X   : A data matrix of observations in which rows represent observations, and columns represent variables.
%          X is used to obtain an initial PCA model using pcaModel.
%   Y   : A data matrix used to validate parameter estimates.
%   k   : Currently used to predetermine the number of components - change
%   later
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

%This source incidates we need to explain most of the variance in Y, not X
%http://white.stanford.edu/~knk/Psych216A/FinalProjects/Grace/html/Psych216aFinalProject.html

%The pls package in R has features to return the variance from X only,
%making me think this is the important one. 

%Cross-validation would be ideal, but may not be possible for RPLS





%A=number of PLS latent vectors to be computed

% beta=PLS regression coef?cients

model.aTX = eps;
model.aEX = eps;
model.aTY = eps;  %for now
model.aEY = eps;  %for now
model.eta = 0.9999; 
if nargin>5,
    varargin=reshape(varargin,2,[]);
    for c=1:size(varargin,2),
        switch varargin{1,c}
            case 'eta'
                model.eta = varargin{2,c};
            case 'aTX'
                model.aTX=varargin{2,c};
            case 'aEX'
                model.aEX=varargin{2,c};
            case 'aTY'
                model.aTY=varargin{2,c};
            case 'aEY'
                model.aEY=varargin{2,c};
           
        end
    end
end

convcrit=1*exp(1)^-10; % convergence criteria for the NIPALS algorithm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PART I
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    model.aG = aG;
    %model.eta = eta;
    model.type = 'RPLS'; 

    n = size(Xt,1);
    model.n = n;
    model.N = n;
    %standardize the data
    model.locX = mean(Xt);
    model.locY = mean(Yt);
    model.Rx = corr(Xt);
    model.Sx = cov(Xt);
    model.Sy = cov(Yt);
    model.kSelectType = 'covmat';
    model.kSelect = 'kSelectCpv';
    stand.Xt = zscore(Xt);
    stand.Yt = zscore(Yt);
    
    %Fit the initial PCA model and set control limits.
    model_0 = plsModel(stand.Xt, stand.Yt, k);%, model.kSelect, model.kSelectType );
    model.P = model_0.P;
    model.k = model_0.k;
    model.W = model_0.W;
    model.WT = model_0.WT; %CHECK THIS OUT R is already correlation % I have changed this to WT
    model.Q = model_0.Q;
    model.LX = model_0.LX;
    model.LY = model_0.LY;
    model.XX = stand.Xt'*stand.Xt;
    model.XY = stand.Xt'*stand.Yt;
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PART II
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if isempty(model.eta) == 1, 
    [model.eta] = fSpeedSolvePls([],Xv, Yv, model);
    end
    %Find optimal alpha values
    [model] = aSolvePls(Xv,Yv, model, aG, model.aTX, model.aEX, model.aTY, model.aEY);
    %Estimate final model using optimal values of alpha
    [~, ~, model] = rplsMonitor(Xv, Yv, model);
end