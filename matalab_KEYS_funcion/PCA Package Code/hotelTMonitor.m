function [ stat, ucl ] = hotelTMonitor( X, model )
%HOTELTMONITOR applies a HOtelling's T^2 model for monitoring the observation of X.
%
% References: 

% Required input arguments: 
%   X     : A data matrix of observations where rows represent observations, and columns represent variables. 
%   model : a structure witht the PCA model parameters.
%
% I/O: [ stat, ucl ] = hotelTMonitor( X, model );
%
% The outputs of HOTELTMONITOR are structures containing:
%   stat.sT   : score statistic.
%   stat.sE   : residual statistic.
%   ucl.uclT  : upper control limit of the score statistic.
%   ucl.uclE  : upper control limit of the residual statistic.
%
% see also pcaModel.

% Monitoring --------------------------------------------------------------

n=size(X,1);

weightT=model.weightT;
stat.sT=zeros(n,1);
for i=1:n,
    stat.sT(i)=X(i,:)*weightT*X(i,:)';% T^2
end

ucl.uclT=ones(n,1)*model.uclT;


%--------------------------------------------------------------------------

end

