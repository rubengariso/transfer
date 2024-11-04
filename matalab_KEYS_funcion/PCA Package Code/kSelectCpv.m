function [ k ] = kSelectCpv( L, cpvThresh )
% KSELECTCPV selects the number of latent variables to be used on PCA 
% models base on the cumulative percentage of variance (CPV) criterion.
% The required percentage of variance explained is set to 95% by default.
%
% References: 
% 
%    Jolliffe, I. T. (2002),
%    Principal Component Analysis. New York, Springer.
%
%    Krzanowski, W. J. and P. Kline (1995),
%    "Cross-Validation for Choosing the Number of Important Components in Principal Component Analysis.",
%    Multivariate Behavioral Research, 30:2, 149-165.
%
% Required input argument:
%   L  : A atrix of singular values corresponding to the covariance matrix of a data sample. 
%
% I/O: [ k ] = kSelectCpv( X )
%
% The output of KSELECTCPV is:
%   k : The selected number of latent variables


vecL=diag(L);
CPV=cumsum(vecL)/sum(vecL);
k=find(CPV>=cpvThresh,1,'first');


end

