function [stat, ucl, model] = rplsMonitor(X,Y,model)

    [n, mX] = size(X);
    stat.sEX = NaN(n,1);
    stat.sTX =  NaN(n,1);
    stat.sEY = NaN(n,1);
    stat.sTY =  NaN(n,1);    
    ucl.sEX =  NaN(n,1);
    ucl.sTX =  NaN(n,1);
    ucl.sEY =  NaN(n,1);
    ucl.sTY =  NaN(n,1);   
    
    for i = 1:n,
    eta_n = model.n/(model.n+1)*model.eta; % weight parameter    
    xi = X(i,:);
    yi = Y(i,:);
    xi_center = xi - model.locX;% mean centering of the new signal
    xi_scale = (xi_center)./sqrt(diag(model.Sx))'; 
    yi_center = yi - model.locY;% mean centering of the new signal
    yi_scale = (yi_center)./sqrt(diag(model.Sy))'; 
   
    %model.T =  xi_scale*model.P(:,1:model.k); % PCA score of the new signal
    %X statistics
    sEX = sum(( xi_scale -  xi_scale*model.P(:,1:model.k)*model.P(:,1:model.k)')'.^2); % Q statistic, in LV, the index i for res is not needed, we will show res in a waveform chart
    sTX = sum(((xi_scale* model.P(:,1: model.k)*inv(model.LX(1: model.k,1: model.k).^.5)).^2),2);
    stat.sEX(i) = sEX;
    stat.sTX(i) = sTX;
    %Y statistics
%     sEY = sum(( yi_scale -  yi_scale*model.Q(:,1:model.k)*model.Q(:,1:model.k)')'.^2); % Q statistic, in LV, the index i for res is not needed, we will show res in a waveform chart
%     sTY = sum(((yi_scale* model.Q(:,1: model.k)*inv(model.LY.^.5)).^2),2);
%     stat.sEY(i) = sEY;
%     stat.sTY(i) = sTY;
    
    %RPLS update begins
    %if sEX<=model.ucl.sEX & sTX<=model.ucl.sTX, % & sEY<=model.ucl.sEY & sTY<=model.ucl.sTY, 
    if i <=n,
        model.locX_old = model.locX;
        model.locX = eta_n*model.locX + (1-eta_n)*xi;
        model.locY_old = model.locY;
        model.locY = eta_n*model.locY + (1-eta_n)*yi;
        
        d_locX = model.locX - model.locX_old; % difference between new and old b
        model.Sx = eta_n*(model.Sx+d_locX'*d_locX)+(1-eta_n)*xi_center'*xi_center; % update for the covariance matrix to be applied for scaling
        d_locY = model.locY - model.locY_old; % difference between new and old b
        model.Sy = eta_n*(model.Sy+d_locY'*d_locY)+(1-eta_n)*yi_center'*yi_center; % update for the covariance matrix to be applied for scaling

        
        if i == 1
            XX = eta_n*model.XX+xi_scale'*xi_scale;
            XY = eta_n*model.XY+xi_scale'*yi_scale;
        else
            XX = eta_n*XX+xi_scale'*xi_scale;
            XY = eta_n*XY+xi_scale'*yi_scale;
        end 
%         xy=xi_scale'*yi_scale; % compute the covariance
%         xx=xi_scale'*xi_scale; % matrices
        M = size(yi, 2);
        M = model.k;
        if i >1 
            clearvars P W Wt Q;
        end
        P = [];
        W = [];
        WT = [];
        Q = [];
            
        for j=1:mX, % A=number of PLS components to be computed
            if M==1, % if there is a single response variable, compute the
                w=XY; % X-weights as shown here
                else % else
                [C,D]=eig(XY'*XY); % ?rst compute the eigenvectors of YTXXTX
                q=C(:,find(diag(D)==max(diag(D)))); % ?nd the eigenvector corresponding to the largest eigenvalue
                w=(XY*q); % compute X-weights
            end
            w=w/sqrt(w'*w); % normalize w to unity
            wt=w; % loop to compute ri
            for h=1:j-1,
                wt=wt-(P(:,h)'*w)*WT(:,h);
            end
            tt=(wt'*XX*wt); % compute tTt
            p=(wt'*XX)'/tt; % X-loadings
            q=(wt'*XY)'/tt; % Y-loadings
            XY=XY-(p*q')*tt; % XTY de?ation
            P=[P p];
            Q=[Q q];
            WT=[WT wt];
            W=[W w];
        end
            model.W=[W w]; % storing loadings and weights
            model.P=P; %original paper code this is [P P]
            model.Q=Q;
            model.WT=WT;
            model.p = p;
            model.q = q;
            model.tt = tt;
            model.LX = diag(var(model.P));
            model.LY = diag(var(model.Q));
            model.ucl.sTX = limitT(model.aTX,model.k);
            model.ucl.sEX = limitQ(model.LX, model.k, model.aEX);
            model.ucl.sTY = limitT(model.aTY,model.k);
            model.ucl.sEY = limitQ(model.LY, model.k, model.aEY);     
            model.n = model.n+1;
    end
        ucl.sEX(i) = model.ucl.sEX;
        ucl.sTX(i) = model.ucl.sTX;
%         ucl.sEY(i) = model.ucl.sEY;
%         ucl.sTY(i) = model.ucl.sTY;    
%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%
    end
    model.beta=model.WT*model.Q'; % compute the regression coefficients
end
