function [ model] = aSolvePls( X,Y, model, aG, aTX, aEX, aTY, aEY )
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
% I/O: [ model] = aSolve( X, model, aG, NaN, NaN );

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

    n = size(X,1);
    funpar=[];
    upper = aG;
    lower = eps;

    if isnan(aTX),
        [ aTX ] = fzero( @(x)funTX(x), [lower upper]);
    end
    if isnan(aEX),
        [ aEX ] = fzero( @(x)funEX(x), [lower upper]);
    end
    if isnan(aTY),
        [ aTY ] = fzero( @(x)funTY(x), [lower upper]);
    end
    if isnan(aEY),
        [ aEY ] = fzero( @(x)funEY(x), [lower upper]);
    end
    
    model.aTX = aTX;
    model.aEX = aEX;
    model.aTY = aTY;
    model.aEY = aEY;
    model.ucl.sTX = limitT(model.aTX,model.k);
    model.ucl.sEX = limitQ(model.LX, model.k, model.aEX);
    model.ucl.sTY = limitT(model.aTY,model.k);
    model.ucl.sEY = limitQ(model.LY, model.k, model.aEY);
           
    function f = funTX(aTX)
        model.aTX = aTX;
        model.aEX = 0;
        model.aTY = 0;
        model.aEY = 0;
        model.ucl.sTX = limitT(model.aTX,model.k);
        model.ucl.sEX = limitQ(model.LX, model.k, model.aEX);
        model.ucl.sTY = limitT(model.aTY,model.k);
        model.ucl.sEY = limitQ(model.LY, model.k, model.aEY);
        
        if strcmp(model.type,'RPLS')
            [stat, ucl] = rplsMonitor(X,Y, model);
        elseif strcmp(model.type,'MWPLS')
            [stat, ucl] = mwplsMonitor(X,Y, model);
        end
        fdrTX = sum(stat.sTX>ucl.sTX)/n;
        f =fdrTX-aG/4;
        
        funpar.fdrTX=fdrTX;

    end

    function f = funEX(aEX)

        model.aTX = 0;
        model.aEX = aEX;
        model.aTY = 0;
        model.aEY = 0;
        model.ucl.sTX = limitT(model.aTX,model.k);
        model.ucl.sEX = limitQ(model.LX, model.k, model.aEX);
        model.ucl.sTY = limitT(model.aTY,model.k);
        model.ucl.sEY = limitQ(model.LY, model.k, model.aEY);        
        
        if strcmp(model.type,'RPLS')
            [stat, ucl] = rplsMonitor(X,Y, model);
        elseif strcmp(model.type,'MWPLS')
            [stat, ucl] = mwplsMonitor(X,Y, model);
        end
        fdrEX = sum(stat.sEX>ucl.sEX)/n;
        f =fdrEX-aG/4;
        
        funpar.fdrEX=fdrEX;
    end

    function f = funTY(aTY)
        model.aTX = 0;
        model.aEX = 0;
        model.aTY = aTY;
        model.aEY = 0;
        model.ucl.sTX = limitT(model.aTX,model.k);
        model.ucl.sEX = limitQ(model.LX, model.k, model.aEX);
        model.ucl.sTY = limitT(model.aTY,model.k);
        model.ucl.sEY = limitQ(model.LY, model.k, model.aEY);
        
        if strcmp(model.type,'RPLS')
            [stat, ucl] = rplsMonitor(X,Y, model);
        elseif strcmp(model.type,'MWPLS')
            [stat, ucl] = mwplsMonitor(X,Y, model);
        end
        fdrTY = sum(stat.sTY>ucl.sTY)/n;
        f =fdrTY-aG/4;
        
        funpar.fdrTX=fdrTX;

    end

    function f = funEY(aEY)

        model.aTX = 0;
        model.aEX = 0;
        model.aTY = 0;
        model.aEY = aEY;
        model.ucl.sTX = limitT(model.aTX,model.k);
        model.ucl.sEX = limitQ(model.L, model.k, model.aEX);
        model.ucl.sTY = limitT(model.aTY,model.k);
        model.ucl.sEY = limitQ(model.L, model.k, model.aEY);        
        
        if strcmp(model.type,'RPLS')
            [stat, ucl] = rplsMonitor(X,Y, model);
        elseif strcmp(model.type,'MWPLS')
            [stat, ucl] = mwplsMonitor(X,Y, model);
        end
        fdrEY = sum(stat.sEY>ucl.sEY)/n;
        f =fdrEY-aG/4;
        
        funpar.fdrEY=fdrEY;
    end


end