function [ model ] = dpcaModel( X, kSelect, varargin )

%DPCAMODEL trains a DPCA model for constructing DPCA control charts. The DPCA
%model performs PCA on an augmented data matrix with time shifted variables 
%selected by lagSelect.
%
% References: 
%
%    Ku, W., R. H. Storer, et al. (1995),
%    "Disturbance detection and isolation by dynamic principal component analysis.",
%    Chemometrics and Intelligent Laboratory Systems 30:1, 179-196.
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
% I/O: model=dpcaModel( X, kSelect, 'maxLag', 10, 'lag', [], 'mode', 'auto' );
%
% The output of DPCAMODEL is a structure containing:
%   model.k       : The selected number of latent variables
%   model.P       : The DPCA model loadings.
%   model.L       : The singular values.
%   model.type    : The type of model being computed is stored as 'DPCA'
%   model.weightT : Weight term in the calculation of the score statistic
%                   (T) used for all evaluations. Stored to save
%                   computation.
%   model.weightE : Weight term in the calculation of the model error
%                   statistic (E), used for all evaluations. Stored to save
%                   computation.
%   model.Xpast   : Past obsevation needed to constrcute the extended
%                   observations vector.
%   model.lag     : The number of lags used for each variable.
%
% see also pcaModel, dpcaMonitor, lagSelect.

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

% extended data matix
Xlag=lagData(X, lag );  

%--------------------------------------------------------------------------

% DPCA model --------------------------------------------------------------

[model]=pcaModel(Xlag,kSelect,'datamat','cpvThresh',cpvThresh,'crossValMode',crossValMode);
model.type='DPCA';
model.lag=lag;

P=model.P(:,1:model.k);
L=model.L(1:model.k,1:model.k);
iL=inv(L);
m=size(Xlag,2);

weightT=P*iL*P';
weightE=eye(m)-P*P';

model.weightT=weightT;
model.weightE=weightE;

model.Xpast=NaN(max(lag),size(X,2));%X(end-max(lag)+1:end,:);

%--------------------------------------------------------------------------


end

