function [ contrT, contrE ] = contributionPlot( input, model )

%PCAMODEL trains a PCA model for constructing PCA control charts. The PCA
%model may be used to monitor in its own right, or as the basis for other
%PCA monitoring methods.
%
% References: 
%    Miller, P., Swanson, R. E., and Heckler, C. H. E. (1998),
%    “Contribution Plots: A Missing Link in Multivariate Quality Control,”
%    Applied Mathematics and Computer Science, 8, 775–792.
%
% Required input arguments: 
%   input   : The observation during which the fault occured whose variables' contributions should be analyzed.  
%   model   : The PCA model used to evaluate the observation at the time of the fault.

%
% I/O: [ contrT, contrE ] = contributionPlot( input, model );
%
% The output of PCAMODEL is a structure containing:
%   contrT        : The selected number of latent variables
%   contrE       : The PCA model loadings.

% see also pcaModel, pcaMonitor

%compute the contributions of each variable
m = length(input);
k = model.k;
iL=inv(model.L(1:k,1:k));
P = model.P(:,1:k);
contrT = input*P*sqrt(iL)*P';
contrE = input*(eye(m)-P*P');

%plot the results
subplot(2,1,1);
bar(contrT)
title('Contribution to T^2-statistic');
subplot(2,1,2);
hold on
bar(contrE)
title('Contribution to Q-statistic');
hold off;
end