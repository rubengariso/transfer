function [ stat, ucl ] = pcaMonitor( X, model )
%PCAMONITOR applies a PCA model for monitoring the observation of X.
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
%   X     : A data matrix of observations where rows represent observations, and columns represent variables. 
%   model : a structure witht the PCA model parameters.
%
% I/O: [ stat, ucl ] = pcaMonitor( X, model );
%
% The outputs of PCAMONITOR are structures containing:
%   stat.sT   : score statistic.
%   stat.sE   : residual statistic.
%   ucl.sT    : upper control limit of the score statistic.
%   ucl.sE    : upper control limit of the residual statistic.
%
% see also pcaModel.

% Monitoring --------------------------------------------------------------

n=size(X,1);

weightT=model.weightT;
weightE=model.weightE;

stat.sT=zeros(n,1);
stat.sE=zeros(n,1);
for i=1:n,
    stat.sT(i)=X(i,:)*weightT*X(i,:)';% T^2
    stat.sE(i)=X(i,:)*weightE*X(i,:)';% Q (SPE)
end

ucl.sT=ones(n,1)*model.ucl.sT;
ucl.sE=ones(n,1)*model.ucl.sE;

%--------------------------------------------------------------------------

end

