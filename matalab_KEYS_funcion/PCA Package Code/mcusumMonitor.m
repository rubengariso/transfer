function [ stat, ucl, model ] = mcusumMonitor( X, model )
% mcusumMonitor
%
% Input arguments:
%     X : Data matrix (observations x variables)
% model : MCUSUM model parameters
%
% Output arguments:
%  stat  : monitoring statistic
%  ucl   : upper control limit
%  model : updated model
%
% I/O: [ stat, ucl, model ] = mcusumMonitor( X, model )

[n,m]=size(X);

stat.mcusum=zeros(n,1);
ucl.mcusum=model.ucl.mcusum*ones(n,1);
for i=1:n,
    d=sqrt((model.C+X(i,:)-model.mew)*model.invS0*(model.C+X(i,:)-model.mew)');
    if d<=model.k,
        model.C=zeros(1,m);
    else
        model.C=(model.C+X(i,:)-model.mew)*(1-model.k/d);
    end
    
    stat.mcusum(i)=sqrt(model.C*model.invS0*model.C');
end


end

