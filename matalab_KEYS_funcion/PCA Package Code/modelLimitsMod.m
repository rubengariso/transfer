function [ model ] = modelLimitsMod( X, model, aG, varargin )
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

Y = NaN; 
if nargin>3,
    varargin=reshape(varargin,2,[]);
    for c=1:size(varargin,2),
        switch varargin{1,c}
            case 'PLS'
                Y = varargin{2,c};
         end
    end
end


if strcmp(model.type,'HotellingT')~=1
    model.uclE=NaN;
end
model.uclT=NaN;

if strcmp(model.type,'PLS')==1
    model.uclEX=NaN;
    model.uclEY=NaN;
    model.uclTX=NaN;
    model.uclTY=NaN;
    
end

if strcmp(model.type,'PCA')==1
    [ stat ] = pcaMonitor( X, model );
elseif strcmp(model.type,'DPCA')==1
    [ stat ] = dpcaMonitor( X, model );
elseif strcmp(model.type,'DPCADR')==1
    [ stat ] = dpcadrMonitor( X, model );
elseif strcmp(model.type,'HotellingT')==1
    [ stat ] = hotelTMonitor( X, model );
elseif strcmp(model.type,'PLS')==1
    [ stat ] = plsMonitor( X,Y, model );
end
if strcmp(model.type,'HotellingT')==1
model.uclT = prctile(stat.sT,(1-aG)*100);
else
    %Set limits for static PCA methods
    if strcmp(model.type,'PLS')~=1
    [ UCLaux, alpha_ind, aGout ] = setLimits( [ stat.sT stat.sE ], aG );
    model.uclT=UCLaux(1);
    model.uclE=UCLaux(2);
    
    %Set limits for static PLS methods
    else
    [ UCLaux, alpha_ind, aGout ] = setLimits( [ stat.sTX stat.sEX ], aG/2 );
    model.uclTX=UCLaux(1);
    model.uclEX=UCLaux(2);
    [ UCLaux, alpha_ind, aGout ] = setLimits( [ stat.sTY stat.sEY ], aG/2 );
    model.uclTY=UCLaux(1);
    model.uclEY=UCLaux(2);
    end
end

end

