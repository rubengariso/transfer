function [ stat, ucl, model ] = dpcadrMonitor( X, model )
%DPCADRMONITOR applies a DPCA model for monitoring the observation of X.
%
% Reference: 
%
%    Rato, T. J. and M. S. Reis (2013),
%    "Fault detection in the Tennessee Eastman benchmark process using dynamic principal components analysis based on decorrelated residuals (DPCA-DR).",
%    Chemometrics and Intelligent Laboratory Systems 125:15, 101-108.
%
% Required input arguments: 
%   X       : A data matrix of observations where rows represent observations, and columns represent variables.  
%   model   : A structure witht the DPCADR model parameters.
%
% I/O: [ stat, ucl, model ] = dpcadrMonitor( X, model );
%
% The outputs of DPCADRMONITOR are structures containing:
%   stat.sT        : score statistic.
%   stat.sE        : residual statistic.
%   ucl.sT         : upper control limit of the score statistic.
%   ucl.sE         : upper control limit of the residual statistic.
%   model.k        : the selected number of latent variables
%   model.P        : the DPCA model loadings.
%   model.L        : the singular values.
%   model.type     : the type of model being computed is stored as 'DPCADR'
%   model.m        : number of variables.
%   model.loc_prev : sample mean of scores residuals.
%   model.S_prev   : sample covariance of scores residuals.
%   model.iS_prev  : inverse if the sample covariance of scores residuals.  Stored to save
%                    computation.
%   model.loc_res  : sample mean of observations residuals.
%   model.S_res    : sample covariance of observations residuals.
%   model.iS_res   : inverse if the sample covariance of observations residuals.  Stored to save
%                    computation.
%   model.Xpast    : past obsevation needed to constrcute the extended
%                    observations vector.
%   model.lag      : the number of lags used for each variable.
%
% see also dpcadrModel.

% pre-processement --------------------------------------------------------

Z=[model.Xpast; X];
n_ext=size(Z,1);
n_lag=max(model.lag);
n_obs=n_ext-n_lag;

% extended data matrix
Zlag = lagData(Z, model.lag );

% correct matrix dimention
Zobs=Z(end-n_obs+1:end,:);
Zobslag=Zlag(end-n_obs+1:end,:);

%--------------------------------------------------------------------------

% DPCA-DR -----------------------------------------------------------------

[ Test ] = missingData( Zobslag, model );
T=Zobslag*model.P(:,1:model.k);
Et=T-Test;
Z_prev=Et-repmat(model.loc_prev,n_obs,1);

Zest=Test*model.P(:,1:model.k)';
Zest=Zest(:,1:model.m);
Ez=Zobs-Zest;
Z_res=Ez-repmat(model.loc_res,n_obs,1);

% monitoring statistics
stat.sT=nan(n_obs,1);%T2_prev
stat.sE=nan(n_obs,1);%T2_res

for i=1:n_obs,
    stat.sT(i)=Z_prev(i,:)*model.iS_prev*Z_prev(i,:)';% T2 to the scores
    stat.sE(i)=Z_res(i,:)*model.iS_res*Z_res(i,:)';% T2 to the residuals
end

ucl.sT=ones(n_obs,1)*model.ucl.sT;
ucl.sE=ones(n_obs,1)*model.ucl.sE;

model.Xpast=Z(end-n_lag+1:end,:);

%--------------------------------------------------------------------------

end

