function [ T, model ] = missingData( X, model )
%MISSINGDATA applys missing data estimation to the first m variables in X
% based on a PCA model.
%
% Reference: 
%
%   Nelson, P. R. C., P. A. Taylor, et al. (1996),
%   "Missing data methods in PCA and PLS: Score calculations with incomplete observations.",
%   Chemometrics and Intelligent Laboratory Systems 35, 45-65.
%
% Required input arguments: 
%   X     : A data matrix of observations where rows represent observations, and columns represent variables. 
%   model : A structure witht the DPCA model parameters.
%
% I/O: [ T, model ] = missingData( X, model );
%
% The output of MISSINGDATA are:
%   T     : Estimated scores.
%   model : Orignal input model with a model.MDmatrix field correspondent
%           to the missing data model. Stored to save computation.


if isfield(model,'MDmatrix')==0,
    
    Pp=model.P(model.m+1:end,:);
    
    a=size(model.L,1);% K - total number of variables
    B=[eye(model.k) zeros(model.k,a-model.k)];
    MDmatrix=B*model.L*Pp'*inv(Pp*model.L*Pp');
    
    model.MDmatrix=MDmatrix;
      
end

Zp=X(:,model.m+1:end);

T=model.MDmatrix*Zp';
T=T';

end

