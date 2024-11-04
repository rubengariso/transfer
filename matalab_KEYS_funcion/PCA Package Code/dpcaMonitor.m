function [ stat, ucl, model ] = dpcaMonitor( X, model )
%DPCAMONITOR applies a DPCA model for monitoring the observation of X.
%
% References: 
%
%    Ku, W., R. H. Storer, et al. (1995),
%    "Disturbance detection and isolation by dynamic principal component analysis.",
%    Chemometrics and Intelligent Laboratory Systems 30:1, 179-196.
%
% Required input arguments: 
%   X     : A data matrix of observations where rows represent observations, and columns represent variables. 
%   model : a structure witht the DPCA model parameters.
%
% I/O: [ stat, ucl, model ] = dpcaMonitor( X, model );
%
% The outputs of DPCAMONITOR are structures containing:
%   stat.sT       : score statistic.
%   stat.sE       : residual statistic.
%   ucl.sT        : upper control limit of the score statistic.
%   ucl.sE        : upper control limit of the residual statistic.
%   model.k       : the selected number of latent variables
%   model.P       : the DPCA model loadings.
%   model.L       : the singular values.
%   model.type    : the type of model being computed is stored as 'DPCA'
%   model.weightT : weight term in the calculation of the score statistic
%                   (T) used for all evaluations. Stored to save
%                   computation.
%   model.weightE : weight term in the calculation of the model error
%                   statistic (E), used for all evaluations. Stored to save
%                   computation.
%   model.Xpast   : past obsevation needed to constrcute the extended
%                   observations vector.
%   model.lag     : the number of lags used for each variable.     
%
% see also dpcaModel.

% Pre-processing ----------------------------------------------------------

X=[model.Xpast; X];

n_ext=size(X,1);
n_lag=max(model.lag);
n_obs=n_ext-n_lag;

% extended data matrix
Xlag = lagData(X, model.lag );

%--------------------------------------------------------------------------

% Monitoring --------------------------------------------------------------

weightT=model.weightT;
weightE=model.weightE;

stat.sT=zeros(n_obs,1);
stat.sE=zeros(n_obs,1);
for i=1:n_obs,
    stat.sT(i)=Xlag(i,:)*weightT*Xlag(i,:)';% T^2
    stat.sE(i)=Xlag(i,:)*weightE*Xlag(i,:)';% Q (SPE)
end

ucl.sT=ones(n_obs,1)*model.ucl.sT;
ucl.sE=ones(n_obs,1)*model.ucl.sE;

model.Xpast=X(end-n_lag+1:end,:);

%--------------------------------------------------------------------------

end

