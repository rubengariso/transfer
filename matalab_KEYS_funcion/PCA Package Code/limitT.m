function [uclT] = limitT(aT, k)
%LIMITT calculates an upper control limit for the T^2 statistic.
%
%    J. Edward Jackson (1959),
%    "Quality Control Methods for Several Related Variables",
%    Technometrics, 1:4, 359-377
%
% Required input arguments: 
%   aT : The desired false detection rate for the T^2 statistic.
%   k  : The selected number of latent variables.
%
% I/O: uclT = limitT(aT, k);
%
% The output of LIMITT is the value:
%   uclT :  the upper control limit for the T^2 statistic.
%
% see also rpcaMonitor, mwpcaMonitor, aSolve, fSpeedSolve, limitT

uclT = chi2inv((1-aT),k);
end