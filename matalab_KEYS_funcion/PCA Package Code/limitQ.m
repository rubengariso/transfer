function [ uclE ] = limitQ(L, k, aE)
%LIMITQ calculates an upper control limit for the Q-statistic.
%Uses the limits of Box (1954) as default, Jackson and Mudholkar (1979) limits are also
%provided in commented code
%
%    Box G. E. P. (1954),
%    "Some Theorems on Quadratic Forms Applied in the Study of Analysis of Variance Problems, I. Effect of Inequality of Variance in the One-Way Classification",
%    The Annals of Mathematical Statistics, 25:2, 290-302
%
%    J. Edward Jackson and Govind S. Mudholkar (1979),
%    "Control Procedures for Residuals Associated with Principal Component Analysis",
%    Technometrics, 21:3, 341-349
%
% Required input arguments: 
%   L  : The singular values of the model the limit is being calculated for. 
%   aE : The desired false detection rate for the T^2 statistic.
%        uclE0 : An the value 
%   k  : The selected number of latent variables.
%
% I/O: uclE = limitQ(aE, k);
%
% The output of LIMITQ is the value:
%   uclE :  the upper control limit for the Q-statistic.
%
% see also rpcaMonitor, mwpcaMonitor, aSolve, fSpeedSolve, limitT

vL=sort(diag(L),'descend');
theta=zeros(1,3);
for i=1:3,
    theta(i)=sum(vL(k+1:end).^i);
end

% % Box (1954)
g=theta(2)/theta(1);
h=theta(1)^2/theta(2);
uclE=g*chi2inv((1-aE),h);

% Jackson and Mudholkar (1979)
% h0=1-(2*theta(1)*theta(3))/(3*theta(2)^2);
% c=norminv(1-aE);
% uclE=theta(1)*(c*sqrt(2*theta(2)*h0^2)/theta(1)+1+theta(2)*h0*(h0-1)/theta(1)^2)^(1/h0);
% 

end

