function [ stat, ucl, model ] = m1z2Monitor( X, model )
%M1Z2MONITOR Summary of this function goes here
%   Detailed explanation goes here

% Reynolds and Cho (2006) -------------------------------------------------

EZ0=model.EZ0;% start as ones(p,1)
% EA0=model.EA0;% start as ones(p,1)
lambda=model.lambda;
iSsq=model.iSsq;% Ssq=Sigma.^2; iSsqrt=inv(Ssq);
% AnZ=model.AnZ;% iS=inv(Sigma); diS=diag(diag(iS)); AnZ=diS^-0.5*iS; 
% iSAsq=par_MEWMA.iSAsq;% iS=inv(Sigma); diS=diag(diag(iS)); SA=diS^-0.5*iS*diS^-0.5; SAsq=SA.^2; iSAsq=inv(SAsq); 


[n,m]=size(X);

stat.m1z2=zeros(n,1);
ucl.m1z2=model.ucl.m1z2*ones(n,1);
% stat.M2Z2=zeros(n,1);
% stat.M1A2=zeros(n,1);
% stat.M2A2=zeros(n,1);
for i=1:n,
    
%     z=AnZ*X(i,:)';
    
    EZ=zeros(m,1);
%     EA=zeros(m,1);
    for j=1:m,
        EZ(j)=(1-lambda)*EZ0(j)+lambda*X(i,j)^2;
%         EA(j)=(1-lambda)*EA0(j)+lambda*z(j)^2;
    end
    
    stat.m1z2(i)=(2-lambda)/(2*lambda)*(EZ-1)'*iSsq*(EZ-1);
%     stat.M2Z2(i)=(2-lambda)/(2*lambda)*EZ'*iSsq*EZ;
%     
%     stat.M1A2(i)=(2-lambda)/(2*lambda)*(EA-1)'*iSAsq*(EA-1);
%     stat.M2A2(i)=(2-lambda)/(2*lambda)*EA'*iSAsq*EA;
    
    % update values
    EZ0=EZ;
% 	EA0=EA;

end

% update model
model.EZ0=EZ;
% model.EA0=EA;

%--------------------------------------------------------------------------

end

