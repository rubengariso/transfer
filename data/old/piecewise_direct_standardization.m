function [F,RMSE,Yhat] = piecewise_direct_standardization(X,Y,windoh)

i_max=size(Y,2);
j_limi=size(X,2);
scaling_method='mean-centering';
F=zeros(j_limi,i_max);

for i=1:i_max
    j_min=i-windoh;
    j_max=i+windoh;
    if j_min<1
       j_min=1;
    end
    if j_max>j_limi
       j_max=j_limi;
    end
    [ kpls ] = kSelectPLS_MCCV_3d(X(:,j_min:j_max),Y(:,i),[],scaling_method );
    meanY=mean(Y);
    stdY=ones(1,size(Y,2));
    zY=zScale(Y,meanY,stdY);
    [ ~,~,~,~,BETA, ~ ] = plsModel(X(:,j_min:j_max), zY(:,i) ,kpls );
    F(j_min:j_max,i)=BETA.Y;
end

Yhat=X*F;
E=Y-Yhat;
RMSE=sqrt(mean(E.^2,2));

end

