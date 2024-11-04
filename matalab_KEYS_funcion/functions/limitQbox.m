function [ucl] = limitQbox(stat,alpha)
% 09.05.2023

mew=mean(stat);
sigma2=var(stat);

g=sigma2/(2*mew);
h=2*mew^2/sigma2;

ucl = g*chi2inv((1-alpha),h);

end

