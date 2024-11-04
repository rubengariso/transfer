function [ model ] = modelLimits( X, model, aG )
% MODELLIMITS determines the control limits for the non-adaptive methods based
% on PCA (PCA, DPCA and DPCADR) by adjustment to the monitoring statistics 
% computed in a validation data set X.
%
% Required input arguments:
%   X     : A data matrix of observations where rows represent observations, and columns represent variables.
%   model : A structure witht the DPCADR model parameters.
%   aG    : Target false detection rate.
%
% I/O: [ UCLp, alfap, alfaout ] = modelLimits( X, model, aG );
%
% The output of MODELLIMITS is the original model structure containing:
%   model.uclT : Upper control limit for the score monitoring statistic.
%   model.uclE : Upper control limit for the score monitoring statistic.
%
% see also pcaModel, dpcaModel, dpcadrModel.

model.ucl.sT=NaN;
model.ucl.sE=NaN;

if strcmp(model.type,'PCA')==1
    [ stat ] = pcaMonitor( X, model );
elseif strcmp(model.type,'DPCA')==1
    [ stat ] = dpcaMonitor( X, model );
elseif strcmp(model.type,'DPCADR')==1
    [ stat ] = dpcadrMonitor( X, model );
end


statnames=fieldnames(stat);
sn=length(statnames);
ST=zeros(size(stat.(statnames{1}),1),sn);
for i=1:sn,
    ST(:,i)=stat.(statnames{i});
end


[ UCLaux, alpha_ind, aGout ] = setLimits( ST, aG );
for i=1:sn,
    model.ucl.(statnames{i})=UCLaux(i);
end

end

