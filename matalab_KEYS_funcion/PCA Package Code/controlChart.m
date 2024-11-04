function [ stat, ucl, model, model_new ] = controlChart( Xtrain, Xtest, modelType, kSelect, aG, varargin )
%CONTROLCHART Summary of this function goes here
%   Detailed explanation goes here

% options -----------------------------------------------------------------

plot_delay=1;
plot_wind=20;
monType='off-line';
lagMode='auto';
maxLag=10;
aT = NaN;
aE = NaN;
h = NaN;
eta = NaN; 
create_model=1;
trainValRatio=50;
cpvThresh = 0.95;
typeLimit = 'Theoretic';

if nargin>5,
    varargin=reshape(varargin,2,[]);
    for c=1:size(varargin,2),
        switch varargin{1,c}
            case 'plot_delay'
                plot_delay=varargin{2,c};
            case 'plot_wind'
                plot_wind=varargin{2,c};
            case 'monType'
                monType=varargin{2,c};
            case 'lagMode'
                lagMode=varargin{2,c};
            case 'model'
                create_model=0;
                model=varargin{2,c};
                modelType=model.type;
            case 'maxLag'
                maxLag=varargin{2,c};
            case 'eta'
                eta = varargin{2,c};
            case 'h'
                h=varargin{2,c};
            case 'aT'
                aT=varargin{2,c};
            case 'aE'
                aE=varargin{2,c};
            case 'kRef'
                kRef=varargin{2,c};
            case 'trainValRatio'
                trainValRatio=varargin{2,c};
            case 'cpvThresh'
                cpvThresh=varargin{2,c};
            case 'typeLimit'
                typeLimit=varargin{2,c};
        end
    end
end

%--------------------------------------------------------------------------

n=size(Xtest,1);

switch typeLimit
    case 'Theoretic'
        
        Xt1=Xtrain;
        Xt2=[];
        
    case 'Empiric'
        
        midPoint=round(size(Xtrain,1)*trainValRatio/100);
        
        Xt1=Xtrain(1:midPoint,:);
        Xt2=Xtrain(midPoint+1:end,:);
        
end

% model and limits --------------------------------------------------------

if create_model==1,
    switch modelType
        case 'PCA'
            model=pcaModel(Xt1, kSelect, 'datamat','cpvThresh',cpvThresh);
            switch typeLimit
                case 'Theoretic'
                    [ uclT] = limitT(aG/2, model.k);
                    [ uclE ] = limitQ(model.L, model.k, aG/2);
                    model.ucl.sT=uclT;
                    model.ucl.sE=uclE;
                case 'Empiric'
                    model=modelLimits( Xt2, model, aG );
            end
        case 'DPCA'
            model=dpcaModel(Xt1, kSelect, 'mode',lagMode,'maxLag',maxLag,'cpvThresh',cpvThresh);
            switch typeLimit
                case 'Theoretic'
                    [ uclT] = limitT(aG/2, model.k);
                    [ uclE ] = limitQ(model.L, model.k, aG/2);
                    model.ucl.sT=uclT;
                    model.ucl.sE=uclE;
                case 'Empiric'
                    model=modelLimits( Xt2, model, aG );
            end
        case 'DPCADR'
            model=dpcadrModel(Xt1, kSelect, 'mode',lagMode,'maxLag',maxLag,'cpvThresh',cpvThresh);
            switch typeLimit
                case 'Theoretic'
                    [ uclT] = limitT(aG/2, model.k);
                    [ uclE ] = limitT(aG/2, model.m);
                    model.ucl.sT=uclT;
                    model.ucl.sE=uclE;
                case 'Empiric'
                    model=modelLimits( Xt2, model, aG );
            end
        case 'RPCA'
            model = rpcaModel(Xt1, Xt2, aG,'eta',eta,'aT',aT,'aE',aE,'cpvThresh',cpvThresh,'typeLimit',typeLimit);
        case 'MWPCA'
            model = mwpcaModel(Xt1, Xt2, aG,'h',h,'aT',aT,'aE',aE,'cpvThresh',cpvThresh,'typeLimit',typeLimit);
        case 'MCUSUM'
            [ model ] = mcusumModel( mean(Xt1), cov(Xt1), kRef );
            model.ucl.mcusum=NaN;
            [ stat ] = mcusumMonitor( Xt2, model );
            [ model.ucl.mcusum ] = setLimits( stat.mcusum, aG );
        case 'MEWMA'
            [ model ] = mewmaModel( eta, cov(Xt1) );
            model.ucl.mewma=NaN;
            [ stat ] = mewmaMonitor( Xt2, model );
            [ model.ucl.mewma ] = setLimits( stat.mewma, aG );
        case 'M1Z2'
            [ model ] = m1z2Model( Xt1, eta );
             model.ucl.m1z2=NaN;
            [ stat ] = m1z2Monitor( Xt2, model );
            [ model.ucl.m1z2 ] = setLimits( stat.m1z2, aG );
        case 'MAXD'
            [ model ] = maxDModel( Xt1, eta );
            model.ucl.maxD=NaN;
            [ stat ] = maxDMonitor( Xt2, model );
            [ model.ucl.maxD ] = setLimits( stat.maxD, aG );
        case 'LRT'
            [ model ] = lrtModel( Xt1, eta );
            model.ucl.lrt=NaN;
            [ stat ] = lrtMonitor( Xt2, model );
            [ model.ucl.lrt ] = setLimits( stat.lrt, aG );
    end
end

%--------------------------------------------------------------------------

% monitor off-line --------------------------------------------------------

if strcmp(monType,'off-line')==1,
    switch modelType
        case 'PCA'
            [ stat, ucl ] = pcaMonitor( Xtest, model );
            model_new=model;
        case 'DPCA'
            [ stat, ucl, model_new ] = dpcaMonitor( Xtest, model );
        case 'DPCADR'
            [ stat, ucl, model_new ] = dpcadrMonitor( Xtest, model );
        case 'RPCA'
            [stat, ucl, model_new] = rpcaMonitor(Xtest, model);
        case 'MWPCA'
            [stat, ucl, model_new] = mwpcaMonitor(Xtest, model);
        case 'MCUSUM'
            [ stat, ucl, model_new ] = mcusumMonitor( Xtest, model );
        case 'MEWMA'
            [ stat, ucl, model_new ] = mewmaMonitor( Xtest, model );
        case 'M1Z2'
            [ stat, ucl, model_new ] = m1z2Monitor( Xtest, model );
        case 'MAXD'
            [ stat, ucl, model_new ] = maxDMonitor( Xtest, model );
        case 'LRT'
            [ stat, ucl, model_new ] = lrtMonitor( Xtest, model );
    end
end

%--------------------------------------------------------------------------

% monitor on-line ---------------------------------------------------------

if strcmp(monType,'on-line')==1,
    model_new=model;
    figure('name','Plot stat on-line','NumberTitle','off')
    for i=1:n,
        time=tic;
        x=Xtest(i,:);
        switch modelType
            case 'PCA'
                [ stat_aux, ucl_aux ] = pcaMonitor( x, model_new );
                model_new=model;
            case 'DPCA'
                [ stat_aux, ucl_aux, model_new ] = dpcaMonitor( x, model_new );
            case 'DPCADR'
                [ stat_aux, ucl_aux, model_new ] = dpcadrMonitor( x, model_new );
            case 'RPCA'
                [stat_aux, ucl_aux, model_new] = rpcaMonitor(x, model_new);
            case 'MWPCA'
                [stat_aux, ucl_aux, model_new] = mwpcaMonitor(x, model_new);
            case 'MCUSUM'
                [ stat_aux, ucl_aux, model_new ] = mcusumMonitor( x, model_new );
            case 'MEWMA'
                [ stat_aux, ucl_aux, model_new ] = mewmaMonitor( x, model_new );
            case 'M1Z2'
                [ stat_aux, ucl_aux, model_new ] = m1z2Monitor( x, model_new );
            case 'MAXD'
                [ stat_aux, ucl_aux, model_new ] = maxDMonitor( x, model_new );
            case 'LRT'
                [ stat_aux, ucl_aux, model_new ] = lrtMonitor( x, model_new );
        end
        if i==1,
            statnames=fieldnames(stat_aux);
            sn=length(statnames);
            stat.(statnames{1})=nan(n,1);
            ucl.(statnames{1})=nan(n,1);
            hFig=cell(sn,1);
        end
        
        time=toc(time);
        for j=1:sn,
            stat.(statnames{j})(i)=stat_aux.(statnames{j});
            ucl.(statnames{j})(i)=ucl_aux.(statnames{j});
            subplot(sn,1,j)
            [ hFig{j} ] = plotChart( hFig{j}, stat.(statnames{j})(i), ucl.(statnames{j})(i), max(plot_delay-time,0), plot_wind );
            ylabel((statnames{j}))
            time=inf;
        end
    end
end

%--------------------------------------------------------------------------

end

