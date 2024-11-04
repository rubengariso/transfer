function [ model ] = m1z2Model( X, lambda )
%M1Z2MODEL Summary of this function goes here
%   Detailed explanation goes here

m=size(X,2);
Sigma=cov(X);

model.type='M1Z2';

model.EZ0=ones(m,1);% inicia em ones(p,1)
% model.EA0=ones(m,1);% inicia em ones(p,1)
model.lambda=lambda;

Ssq=Sigma.^2;
iSsq=inv(Ssq);
model.iSsq=iSsq;% Ssq=Sigma.^2; iSsq=inv(Ssq);

% iS=inv(Sigma);
% diS=diag(diag(iS));
% AnZ=diS^-0.5*iS;
% model.AnZ=AnZ;% iS=inv(Sigma); diS=diag(diag(iS)); AnZ=diS^-0.5*iS;
% SA=diS^-0.5*iS*diS^-0.5;
% SAsq=SA.^2;
% iSAsq=inv(SAsq);
% model.iSAsq=iSAsq;% iS=inv(Sigma); diS=diag(diag(iS)); SA=diS^-0.5*iS*diS^-0.5; SAsq=SA.^2; iSAsq=inv(SAsq); 

end

