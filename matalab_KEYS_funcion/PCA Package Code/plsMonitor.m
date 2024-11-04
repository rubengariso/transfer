function [ stat, ucl ] = plsMonitor( X,Y, model )
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
%   ucl.uclT  : upper control limit of the score statistic.
%   ucl.uclE  : upper control limit of the residual statistic.
%
% see also pcaModel.

% Monitoring --------------------------------------------------------------

n=size(X,1);

weightTX=model.weightTX;
weightTY=model.weightTY;
weightEX=model.weightEX;
weightEY=model.weightEY;

stat.sTX=zeros(n,1);
stat.sTY=zeros(n,1);
stat.sEX=zeros(n,1);
stat.sEY=zeros(n,1);

for i=1:n,
    stat.sTX(i)=X(i,:)*weightTX*X(i,:)';% T^2 X
    stat.sTY(i)=Y(i,:)*weightTY*Y(i,:)';% T^2 Y
    stat.sEX(i)=X(i,:)*weightEX*X(i,:)';% Q (SPE) X
    stat.sEY(i)=Y(i,:)*weightEY*Y(i,:)';% Q (SPE) Y
end

ucl.uclTX=ones(n,1)*model.uclTX;
ucl.uclTY=ones(n,1)*model.uclTY;
ucl.uclEX=ones(n,1)*model.uclEX;
ucl.uclEY=ones(n,1)*model.uclEY;
end

