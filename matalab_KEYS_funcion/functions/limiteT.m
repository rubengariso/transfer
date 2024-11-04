function [UCL] = limiteT(n,k,alpha)

% n- # de observações 
% k= # de variaveis latentes
% alpha confiança

g= ((n-1)*(n+1)*k) / (n* (n - k ));

UCL=g* finv((1-alpha), k, (n-k));
end