function [ Xlag ] = lagData( X, lag )
% LAGDATA creates an extended data matrix by adding lags to each variable
% of X as difined by lag.
%
% Required input arguments:
%   X   : A data matrix of observations where rows represent observations, and columns represent variables. 
%   lag : Vector with the number of lags for each variable in X.
%
% I/O: [ Xlag ] = lagData( X, lag )
%
% The output of LAGDATA is:
%   Xlag : A extended data matrix with lagged variables.

[n, m] =size(X);
d=max(lag);
Xlag=zeros(n-d,m*(d+1));


for i=1:d+1,
    Xlag(:,m*(d+1-i)+1:m*(d+2-i))=X(i:n-d+i-1,:);
end
    
if length(lag)>1,
    % remove unnecessary lags
    in_lag=true(m,d+1);% row: variable, column: lag of each variable
    for i=1:m,
        in_lag(i,lag(i)+2:end)=false;
    end
    Xlag=Xlag(:,in_lag(:));
end

end
