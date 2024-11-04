function [ k ] = kSelectParAnl( X )
% KSELECTPARANL selects the number of latent variables to be used on PCA 
% models base on parallel analysis.
%
% Required input argument:
%   X : A data matrix of observations where rows represent observations, and columns represent variables. 
%
% I/O: [ k ] = kSelectParAnl( X )
%
% The output of KSELECTPARANL is:
%   k : The selected number of latent variables

[n,m]=size(X);
A=random('normal',0,1,n,m);
LX=svd(cov(X));
LA=svd(cov(A));
Ldif=LX-LA;
index=find(Ldif<0,1,'first');
k=index-1;

end

