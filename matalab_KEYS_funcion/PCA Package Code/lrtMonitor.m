function [ stat, ucl, model ] = lrtMonitor( X, model )
%CMONITOR Summary of this function goes here
%   Detailed explanation goes here

[n,m]=size(X);
lambda=model.lambda;
AnU=model.AnU;
Sn0=model.Sn0;
mew_x=model.mew_x;

stat.lrt=zeros(n,1);
ucl.lrt=model.ucl.lrt*ones(n,1);
% Hawkins (2008) ----------------------------------------------------------

for i=1:n,
    
    u=AnU*(X(i,:)-mew_x)';
    Sn=(1-lambda)*Sn0+lambda*(u*u');
    
    stat.lrt(i)=trace(Sn)-log(det(Sn))-m;
    
    % update values
    Sn0=Sn;
end


% update model
model.Sn0=Sn;% inicia em Sn0=eye(p)

%--------------------------------------------------------------------------

end

