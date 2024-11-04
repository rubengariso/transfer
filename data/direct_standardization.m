function [F_DS,model,RMSE,Yhat] = direct_standardization(X_REF,X_NEW,index_mccv)
% X_REF is the domain we want to project to
%X_new is the original domain
Y=X_REF;
X=X_NEW;
scaling_method='mean-centering';
[ kpls ] = kSelectPLS_MCCV_3d(X,Y,index_mccv,scaling_method );

% meanY=mean(Y);
% stdY=ones(1,size(Y,2));
% zY=zScale(Y,meanY,stdY);

model.k=kpls;
[ model.P,  model.Q,  model.B,  model.W,  model.BETA, ~ ] = plsModel(X_NEW ,Y, kpls );
F_DS=model.BETA.Y;

Yhat=X_NEW*F_DS;

E=Y-Yhat;
RMSE=sqrt(mean(E.^2,2));
end

