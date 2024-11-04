function [ model ] = mewmaModel( lambda, S0 )
% mewmaModel
%
% Input arguments:
% lambda : forgeting factor
%     S0 : covariance matrix
%
% Output arguments:
%  model : MEWMA model
%
% I/O: [ model ] = mewmaModel( lambda, S0 )


m=size(S0,2);

model.type='MEWMA';
model.lambda=lambda;
model.S0=S0;
model.invS0=inv(S0);
model.t=0;
model.z=zeros(1,m);

end

