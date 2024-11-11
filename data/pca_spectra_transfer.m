function [Xtransf,RMSE] = pca_spectra_transfer(X_source,X_destination,index_mccv)
%11.11.2024


X=[X_source,X_destination];
m=size(X_source,2);

numMCCV=size(index_mccv.train,2);

RMSE=nan(numMCCV,1);
parfor i=1:numMCCV,

[RMSE(i)]=fun(X(index_mccv.train(:,i),:),X(index_mccv.test(:,i),:),m);
end

[~,Xtransf]=fun(X,X,m);

end

function [RMSE,Xtransf]=fun(Xtrain,Xtest,m)
kSelect='kSelectCpv';
[model]=pcaModel(Xtrain,kSelect,'datamat');
model.m=m;

%--------------------------------------------------------------------------
MDmethod='TSR';
[ T, model ] = missingData( Xtest, model, MDmethod );% estimated scores
Xrec=T*model.P(:,1:model.k)';
Xtransf=Xrec(:,model.m+1:end);
E=Xtest(:,m+1:end)-Xtransf;
RMSE=sqrt(mean(E(:).^2));
end
