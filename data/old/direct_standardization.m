function [F_DS,model,RMSE,Yhat] = direct_standardization_2(X_REF,X_NEW)
% X_REF is the domain we want to project to
%X_new is the original domain
Y=X_REF;
X=X_NEW;
i_total=size(Y,2);

for i=1:i_total,
    [ kpls ] = kSelectPLS_MCCV(X,Y(:,i) );
    
     meanY=mean(Y(:,i));
    stdY=std(Y(:,i));
    zY=zScale(Y(:,i),meanY,stdY);
    
    model.k=kpls;
    [ model.P,  model.Q,  model.B,  model.W,  model.BETA, ~ ] = plsModel(X,zY, kpls );

    F_DS(:,i)=model.BETA.Y;


end

Yhat=X*F_DS;

E=Y-Yhat;
RMSE=sqrt(mean(E.^2,2));
end

