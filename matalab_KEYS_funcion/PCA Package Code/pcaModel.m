function [ model ] = pcaModel( input, kSelect, type, varargin )

%PCAMODEL trains a PCA model for constructing PCA control charts. The PCA
%model may be used to monitor in its own right, or as the basis for other
%PCA monitoring methods.
%
% References: 
%    J. Edward Jackson (1959),
%    "Quality Control Methods for Several Related Variables",
%    Technometrics, 1:4, 359-377
%
%    J. Edward Jackson and Govind S. Mudholkar (1979),
%    "Control Procedures for Residuals Associated with Principal Component Analysis",
%    Technometrics, 21:3, 341-349
%
% Required input arguments: 
%   input   : A covariance matrix or a data matrix of observations. 
%             If Input is a data matrix, then rows of Input represent observations, and columns represent variables.  
%   kSelect : Indicates whether the CPV ('kSelectCpv'), crossvalidation ('kSelectCrossVal') or parallel
%             analysis ('kSelectParAnl') will be used to select the number of latent variables.
%   type    : Indicates whether the input is a covariance matrix ('covmat')or a data
%             matrix ('datamat').
%
% Optional input arguments:
%   cpvThresh    : Threshold for the cumulative percentage of variance used
%                  by CPV ('kSelectCpv')
%   crossValMode : Method used by crossvalidation ('kSelectCrossVal') to determine the PRESS, wheiter 5-fold ('Kfold') or apriximation('aprox').
%
%
% I/O: model=pcaModel( input, kSelect, type );
%
% The output of PCAMODEL is a structure containing:
%   model.k       : The selected number of latent variables
%   model.P       : The PCA model loadings.
%   model.L       : The singular values.
%   model.type    : The type of model being computed is stored as 'PCA'
%   model.weightT : Weight term in the calculation of the score statistic
%                   (T) used for all evaluations. Stored to save
%                   computation.
%   model.weightE : Weight term in the calculation of the model error
%                   statistic (E), used for all evaluations. Stored to save
%                   computation.
%
% see also pcaMonitor


cpvThresh = 0.95;
crossValMode = 'aprox';
if nargin>3,
    varargin=reshape(varargin,2,[]);
    for c=1:size(varargin,2),
        switch varargin{1,c}
            case 'cpvThresh'
                cpvThresh=varargin{2,c};
            case 'crossValMode'
                crossValMode=varargin{2,c};
        end
    end
end

if strcmp(type, 'covmat')
   covmat = input; % singular value decomposition of the covariance matrix, is named SVD factorization in L
elseif strcmp(type, 'datamat')
    X = input;
    covmat = cov(X);
end
[~,model.L,model.P] = svd(covmat);    

if  strcmp(kSelect,'kSelectCpv')==1
    model.k = kSelectCpv( model.L, cpvThresh);
    model.cpvThresh=cpvThresh;
elseif  strcmp(kSelect,'kSelectCrossVal')==1,
    model.k = kSelectCrossVal( X, crossValMode );
elseif strcmp(kSelect,'kSelectParAnl')==1
    model.k = kSelectParAnl( X );
end

% PCA model for monitoring ------------------------------------------------

model.type='PCA';

P=model.P(:,1:model.k);
L=model.L(1:model.k,1:model.k);
iL=inv(L);
m=size(input,2);

weightT=P*iL*P';
weightE=eye(m)-P*P';

model.weightT=weightT;
model.weightE=weightE;

%--------------------------------------------------------------------------
end
