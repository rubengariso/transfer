function [ model] = aSolve( X, model, aG, aT, aE )
%ASOLVE automatically determines the values alpha for the T^2 and Q
%statistic. To find the aT, aE is set to 0 so that no Q-statistic faults
%are detected, and a value of aT is identified that attains the desired
%false detection rate. The reverse approach is used to for finding a
%suitable value of aE.
%
% Required input arguments: 
%   X     : A data matrix of observations in which rows represent observations, and columns represent variables.
%   model : An 'RPCA' or 'MWPCA' model structure to start monitoring from.
%   aG    : The desired global false detection rate.
%   
%
% Optional input arguments:
%   aT : A user selected value for the alpha parameter of the T^2-
%        statistic.
%   aE : A user selected value for the alpha parameter of the Q-
%        statistic.
%
% I/O: [ model] = aSolve( X, model, aG, aT, aE );

%
% The output of ASOLVE is a model structure containing:
%   model :  ASOLVE will return the parameters of the RPCA or MWPCA model
%            that was used as input. Additionally it will provide selected values
%            for the alpha parameters and their limits. These are the
%            following:
%   model.aT : The selected value of aT.
%   model.aE : The selected value of aE.
%   model.ucl.sT : The upper control limit corresponding to the selected value of aT.
%   model.ucl.sE : The upper control limit corresponding to the selected value of aE.
%
% see also rpcaModel, mwpcaModel

    if nargin<=3,
        aT=NaN;
    end
    if nargin<=4,
        aE=NaN;
    end
    
    

    n = size(X,1);
    funpar=[];
    %upper = aG;
    upper = 1-eps;
    lower = eps;

    if isnan(aT),
        [ aT ] = fzero( @(x)funT(x), [lower upper]);
    end
    if isnan(aE),
        [ aE ] = fzero( @(x)funE(x), [lower upper]);
    end
    
    model.aT = aT;
    model.aE = aE;
    model.ucl.sT = limitT(model.aT,model.k);
    model.ucl.sE = limitQ(model.L, model.k, model.aE);
           
    function f = funT(aT)
        model.aT = aT;
        model.aE = 0;
        model.ucl.sT = limitT(model.aT,model.k);
        model.ucl.sE = limitQ(model.L, model.k, model.aE);
        
        if strcmp(model.type,'RPCA')
            [stat, ucl] = rpcaMonitor(X, model);
        elseif strcmp(model.type,'MWPCA')
            [stat, ucl] = mwpcaMonitor(X, model);
        end
        fdrT = sum(stat.sT>ucl.sT)/n;
        f =fdrT-aG/2;
        
        funpar.fdrT=fdrT;

    end

    function f = funE(aE)

        model.aT = 0;
        model.aE = aE;
        model.ucl.sT = limitT(model.aT,model.k);
        model.ucl.sE = limitQ(model.L, model.k, model.aE);
        
        if strcmp(model.type,'RPCA')
            [stat, ucl] = rpcaMonitor(X, model);
        elseif strcmp(model.type,'MWPCA')
            [stat, ucl] = mwpcaMonitor(X, model);
        end
        fdrE = sum(stat.sE>ucl.sE)/n;
        f =fdrE-aG/2;
        
        funpar.fdrE=fdrE;
    end


end