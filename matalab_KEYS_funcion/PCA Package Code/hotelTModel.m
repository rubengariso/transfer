function [ model ] = hotelTModel( input, type )

%HOTELTMODEL trains a Hotelling's T^2 model for constructing Hotelling's T^2 control charts.
%This is performed by fitting a PCA model with 100% of the variation
%explained.
%
% References: 

% Required input arguments: 
%   input   : A covariance matrix or a data matrix of observations. 
%             If Input is a data matrix, then rows of Input represent observations, and columns represent variables.  
%   type    : Indicates whether the input is a covariance matrix ('covmat')or a data
%             matrix ('datamat').
%
% I/O: model=hotelTModel( input, type );
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
% see also hotelTMonitor

if strcmp(type, 'covmat')
   covmat = input; % singular value decomposition of the covariance matrix, is named SVD factorization in L
elseif strcmp(type, 'datamat')
    X = input;
    covmat = cov(X);
end
[~,model.L,model.P] = svd(covmat);    

model.k = size(input, 2);

    % Hotelling's T^2 model for monitoring ------------------------------------------------

model.type='HotellingT';

P=model.P(:,1:model.k);
L=model.L(1:model.k,1:model.k);
iL=inv(L);
m=size(input,2);

weightT=P*iL*P';

model.weightT=weightT;

%--------------------------------------------------------------------------
end
