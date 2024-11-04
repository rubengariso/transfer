function [model] = mwpcaModel(Xt, Xv, aG, varargin)
%MWPCAMODEL outputs an MWPCA model for subsequent online monitoring. This
%function can also perform parameter selection for the the alpha and window
%size parameters. The user can also specify the parameter values instead of
%automatically selecting them.
%
% References: 
%     Wang, X., Kruger, U., Irwin, G. W.,
%     "Process Monitoring Approach Using Fast Moving Window PCA",
%     Industrial & Engineering Chemistry Research, 15, 5691-5702
%
%
% Required input arguments: 
%   Xt   : A data matrix of observations in which rows represent observations, and columns represent variables.
%          Xt is used to obtain an initial PCA model using pcaModel.
%   Xv   : A data matrix used to validate parameter estimates.
%   aG   : The desired global false detection rate.
%
% Optional input arguments:
%   h  : The user provided value of the window size h.
%   aT : The user provided value of the T^2-statistic alpha
%            parameter. aT takes values such that 0 <= aT <= 1.
%   aE : The user provided value of the model Q-statistic alpha
%            parameter. aE takes values such that 0 <= aT <= 1.
%
% I/O: model=mwpcaModel( Xt, Xv, aG, 'h', [], 'aT', [],'aE', []);
%
% The output of MWPCAMODEL is a structure containing:
%   model.type        : The type of model being computed is stored as 'MWPCA'.
%   model.aT          : Value of the T^2-statistic alpha parameter.
%   model.aE          : Value of the Q-statistic alpha parameter.
%   model.h           : Value of h, the window size parameter.
%   model.n           : The number of observations from the dataset Xt used
%                       to train the initial PCA model. Typically, this equals h, but may be
%                       less if a selected h is larger than the number of observations in Xt.
%   model.loc         : The weighted sample mean.
%   model.R           : The weighted sample correlation matrix.
%   model.sd          : The weighted sample standard deviations.
%   model.kSelectType : Latent variable selection always performed using
%                       'covmat' as the input type.
%   model.kSelect     : Latent variable selection always performed using
%                       'CPV' as the seelction method.
%   model.P           : The PCA model loadings.
%   model.L           : The singular values.
%   model.k           : Latent variable selection always performed using 'covmat'
%   model.ucl.sT      : The upper control limit of the T^2-statistic.
%   model.ucl.sE      : The upper control limit of the Q-statistic.
%   model.aG          : The desired global false detection rate
%   model.T           : The current score vector.
%   model.Xh          : The observations in the current, updated window.
%   model.Xh_old      : Used to store observations from the previous period
%                       window.
%
% see also pcaModel, mwpcaMonitor, aSolve, fSpeedSolve.


aT = NaN;
aE = NaN;
h = NaN; 
cpvThresh = 0.95;
typeLimit = 'Theoretic';
if nargin>3,
    varargin=reshape(varargin,2,[]);
    for c=1:size(varargin,2),
        switch varargin{1,c}
            case 'h'
                h = varargin{2,c};
            case 'aT'
                aT=varargin{2,c};
            case 'aE'
                aE=varargin{2,c};
            case 'cpvThresh'
                cpvThresh=varargin{2,c};   
            case 'typeLimit'
                typeLimit=varargin{2,c};
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PART I
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    model.type = 'MWPCA'; 
    model.cpvThresh = cpvThresh;
    if isempty(h) == 1,
        h = fSpeedSolve(Xt,Xv, model, @(Xt, fSpeed, cpvThresh)mwpcaInitialize(Xt, fSpeed, cpvThresh));
    end
    
    [model] = mwpcaInitialize(Xt, h, cpvThresh); 
    model.aG = aG;
    model.cpvThresh = cpvThresh;
    switch typeLimit
        case 'Theoretic'
            modelb.aT = aG/2;
            modelb.aE = aG/2;
        case 'Empiric'
            %Find optimal alpha values
            % [modelb, funpar] = aSolve(Xv, model, aG);
            [modelb] = aSolve(Xv, model, aG, aT, aE);
    end
    model.aT = modelb.aT;
    model.aE = modelb.aE;
    model.ucl.sT = limitT( model.aT,model.k);
    model.ucl.sE = limitQ(model.L, model.k,  model.aE);

    %Estimate final model using optimal values of alpha
    [~, ~, model] = mwpcaMonitor(Xv, model);    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PART II
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    function [model] = mwpcaInitialize(Xt, fSpeed, cpvThresh)
    model.type = 'MWPCA';
    model.aT = 0;
    model.aE = 0;    
    n = size(Xt,1); 
    model.h = round(fSpeed); 
    if n>model.h,
        Xt_use = Xt((n-model.h+1):end,:);
    else
        Xt_use = Xt;
    end
    [n,~] = size( Xt_use);
    model.n = n;
    %standardize the data
    model.loc = mean( Xt_use);
    model.R = corr(Xt_use);
    model.sd = std(Xt_use);
    model.kSelectType = 'covmat';
    model.kSelect = 'kSelectCpv';
    model.cpvThresh = cpvThresh;

    %Fit the initial PCA model and set control limits.
    model_0 = pcaModel(model.R, model.kSelect, model.kSelectType, 'cpvThresh', cpvThresh );
    model.P = model_0.P;
    model.L = model_0.L;
    model.k = model_0.k;
    model.Xh_old = Xt_use; %store the obsevations in the window

    model.ucl.sT = limitT( model.aT, model.k);
    model.ucl.sE = limitQ( model.L, model.k, model.aE);
 
    end

end
