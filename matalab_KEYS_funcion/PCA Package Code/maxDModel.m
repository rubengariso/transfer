function [ model ] = maxDModel( X, lambda )
%MAXDMODEL Summary of this function goes here
%   Detailed explanation goes here

m=size(X,2);

model.type='MAXD';

model.S0z=eye(m);% inicia em S0=eye(p);
model.lambda=lambda;
Sigma=cov(X);
iS=inv(Sigma);
diS=diag(diag(iS));
AnZ=diS^-0.5*iS; 
model.AnZ=AnZ;% iS=inv(Sigma); diS=diag(diag(iS)); AnZ=diS^-0.5*iS; 

model.mewD1=2*m*lambda/(2-lambda);
model.sigmaD1=sqrt(m*(48*lambda^4/(1-(1-lambda)^4)+8*lambda^2/(2-lambda)^2));
model.mewD2=m*(m-1)/2*lambda/(2-lambda);
model.sigmaD2=sqrt(m*(m-1)*(2*m-1)*lambda^4/(1-(1-lambda)^4)+m*(m-1)*(lambda/(2-lambda))^2);


end

