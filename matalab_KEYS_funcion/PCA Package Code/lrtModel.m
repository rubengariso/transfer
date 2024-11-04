function [ model ] = lrtModel( X, lambda )
%CMODEL Summary of this function goes here
%   Detailed explanation goes here

m=size(X,2);
Sigma=cov(X);

model.type='LRT';

model.lambda=lambda;
L=chol(Sigma,'lower');
AnU=inv(L);
model.AnU=AnU;% L=chol(Sigma,'lower'); AnU=inv(L);

model.Sn0=eye(m);% start as Sn0=eye(p)
model.mew_x=mean(X);

end

