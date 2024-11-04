function [ stat, ucl, model ] = mewmaMonitor( X, model )
% mewmaMonitor
%
% Input arguments:
%     X : Data matrix (observations x variables)
% model : MEWMA model parameters
%
% Output arguments:
%  stat : monitoring statistic
%  ucl   : upper control limit
%  model : updated model
%
% I/O: [ stat, ucl, model ] = mewmaMonitor( X, model )

n=size(X,1);

stat.mewma=zeros(n,1);
ucl.mewma=model.ucl.mewma*ones(n,1);
for i=1:n,
    model.t=model.t+1;
    model.z=model.lambda*X(i,:)+(1-model.lambda)*model.z;
    model.c=(model.lambda/(2-model.lambda))*(1-(1-model.lambda)^(2*model.t));

    stat.mewma(i)=model.z*(1/model.c*model.invS0)*model.z';
end



end

