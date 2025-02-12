function [ model ] = dpcadrModel( X, kSelect, varargin )

%DPCADRMODEL trains DPCA and missing data models for constructing DPCADR
% control charts. The DPCADR model performs PCA on an augmented data matrix 
% with time shifted variables selected by lagSelect. The missing data models
% is then trained assuming that the more recent observations are missing.
%
% Reference: 
%
%    Rato, T. J. and M. S. Reis (2013),
%    "Fault detection in the Tennessee Eastman benchmark process using dynamic principal components analysis based on decorrelated residuals (DPCA-DR).",
%    Chemometrics and Intelligent Laboratory Systems 125:15, 101-108.
%
% Required input arguments: 
%   X       : A data matrix of observations where rows represent observations, and columns represent variables.  
%   kSelect : Indicates whether the CPV ('kSelectCpv'), crossvalidation ('kSelectCrossVal') or parallel
%             analysis ('kSelectParAnl') will be used to select the number of latent variables.
%
% Optional input arguments:
%   maxLag       : Maximum number of lag that the lag selection method can use.
%   lag          : Number of lags for each variable.
%   mode         : Lag selection mode, whether automatic ('auto') or manual ('manual').
%   cpvThresh    : Threshold for the cumulative percentage of variance used
%                  by CPV ('kSelectCpv')
%   crossValMode : Method used by crossvalidation ('kSelectCrossVal') to determine the PRESS, wheiter 5-fold ('Kfold') or apriximation('aprox').
%
% I/O: model=dpcadrModel( X, kSelect, 'maxLag', 10, 'lag', [], 'mode', 'auto');
%
% The output of DPCADRMODEL is a structure containing:
%   model.k        : The selected number of latent variables
%   model.P        : The DPCA model loadings.
%   model.L        : The singular values.
%   model.type     : The type of model being computed is stored as 'DPCADR'
%   model.m        : Number of variables.
%   model.loc_prev : Sample mean of scores residuals.
%   model.S_prev   : Sample covariance of scores residuals.
%   model.iS_prev  : Inverse if the sample covariance of scores residuals.  Stored to save
%                    computation.
%   model.loc_res  : Sample mean of observations residuals.
%   model.S_res    : Sample covariance of observations residuals.
%   model.iS_res   : Inverse if the sample covariance of observations residuals.  Stored to save
%                    computation.
%   model.Xpast    : Past obsevation needed to constrcute the extended
%                    observations vector.
%   model.lag      : The number of lags used for each variable.
%   model.MDmatrix : The missing data model. Stored to save computation.
%
% see also pcaModel, dpcadrMonitor, lagSelect, missingData.

% Optional values ---------------------------------------------------------

Lmax=10;
mode='auto';
lag=[];
cpvThresh = 0.95;
crossValMode = 'aprox';

if nargin>2,
    varargin=reshape(varargin,2,[]);
    for c=1:size(varargin,2),
        switch varargin{1,c}
            case 'maxLag'
                Lmax=varargin{2,c};
            case 'lag'
                lag=varargin{2,c};
            case 'mode'
                mode=varargin{2,c};
            case 'cpvThresh'
                cpvThresh=varargin{2,c};
            case 'crossValMode'
                crossValMode=varargin{2,c};
        end
    end
end

%--------------------------------------------------------------------------

% Select number of lags ---------------------------------------------------

if isempty(lag)==1,
    [ lag ] = lagSelect( X, [], Lmax, [], mode );
end

% data dimension
n_ext=size(X,1);
n_lag=max(lag);
n_obs=n_ext-n_lag;
m=size(X,2);

% extended data matix
Xobs=X(end-n_obs+1:end,:);
Xlag=lagData(X, lag );  

%--------------------------------------------------------------------------
    
% DPCA model --------------------------------------------------------------

[model]=pcaModel(Xlag,kSelect,'datamat','cpvThresh',cpvThresh,'crossValMode',crossValMode);
model.type='DPCADR';
model.m=m;

%--------------------------------------------------------------------------

% model for scores
[Test, model]=missingData( Xlag, model );% estimated scores
T=Xlag*model.P(:,1:model.k);% observed scores

E=T-Test;
loc_prev=mean(E);
S_prev=cov(E);
iS_prev=pinv(S_prev);

%--------------------------------------------------------------------------

% model for residuals -----------------------------------------------------

Xest=Test*model.P(:,1:model.k)';
Xest=Xest(:,1:model.m);
E=Xobs-Xest;
loc_res=mean(E);
S_res=cov(E);
iS_res=pinv(S_res);

%--------------------------------------------------------------------------

% store model paramters ---------------------------------------------------

% DPCA model
model.lag=lag;
% T2_prev - scores
model.loc_prev=loc_prev;
model.S_prev=S_prev;
model.iS_prev=iS_prev;
% T2_res - residuals
model.loc_res=loc_res;
model.S_res=S_res;
model.iS_res=iS_res;
model.Xpast=NaN(max(lag),m);%X(end-max(lag)+1:end,:);

%--------------------------------------------------------------------------

end

