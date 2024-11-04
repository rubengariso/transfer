function [ model ] = mcusumModel( mew, S, k )
% mcusumModel
%
% Input arguments:
%   mew : mean of X
%     S : covariance matrix of X
%     k : reference value
%
% Output arguments:
% model : MCUSUM model
%
% I/O: [ model ] = mcusumModel( mew, S, k )


m=size(S,2);

model.type='MCUSUM';
model.mew=mew;
model.S=S;
model.invS0=inv(S);
model.k=k;
model.C=zeros(1,m);

end

