function [F_DS,model,RMSE,Yhat] = direct_standardization(X_source,X_destination,index_mccv)
% X_REF is the domain we want to project to
%X_new is the original domain
% Y=X_REF;
% X=X_NEW;
scaling_method='mean-centering';
[ kpls, RMSE_k, SE_k, n_SE_k   ] = kSelectPLS_MCCV_3d(X_source,X_destination,index_mccv,scaling_method );

% meanY=mean(Y);
% stdY=ones(1,size(Y,2));
% zY=zScale(Y,meanY,stdY);

model.k=kpls;
[ model.P,  model.Q,  model.B,  model.W,  model.BETA, ~ ] = plsModel(X_source ,X_destination, kpls );
F_DS=model.BETA.Y;

Yhat=X_source*F_DS;

% E=Y-Yhat;
% RMSE=sqrt(mean(E.^2,2));
RMSE=RMSE_k(:,end); 
% RMSE_check=sqrt(sum(SE_k(:,end),2)./n_SE_k(:,end)); 
end

