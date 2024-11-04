function [ h, pValue, T0 ] = permtest( y1, y2, test_type, alpha, imax )
% PERMTEST performs a permutation test.
%
% References: 
%     Pesarin, F. and L. Salmaso, Permutation Tests for Complex Data: Theory, 
%         Applications and Software. 2010: John Wiley & Sons Ltd. pp. 19-20.
%
% Required input arguments: 
%   y1     : A column vector of samples.  
%   y2     : A column vector of samples. 
%
% I/O: [ P, Q, B, W ] = plsModel( X, Y, k )
%
% The outputs of PLSMODEL are:
%   h      : result of the hypothesis test, h=0 indicates that the null hypothesis
%            cannot be rejected at the alpha (default 0.05) level. h=1 indicates
%            that the null hypothesis can be rejected at the alpha level.
%   pValue : p-value of the hypothesis test.
%   T0     : test statistic.


% default values ----------------------------------------------------------

if nargin<=2,
    test_type='both';
end

if nargin<=3,
    alpha=0.05;
end

if nargin<=5,
    imax=10000;
end

%--------------------------------------------------------------------------

% permutation test --------------------------------------------------------

% compute diferences and remove missing data
X=y1-y2;
ind=isnan(X)==0;
X=X(ind);
n=length(X);

% generate random permutations
if imax<=10000,
    randPerm=(randi(2,imax,n) > 1)*2 - 1;
end


% compute H0 distribution
T0=sum(X);
if imax<=10000,
    T=randPerm*X;
else
    T=zeros(imax,1);
    for i=1:imax,
        randPerm=(randi(2,1,n) > 1)*2 - 1;
        T(i)=randPerm*X;
    end
end

% performe hypotesis test
switch test_type
    case 'left'
        % p=Pr(T<=t | H0) for a one-sided left-tail test-statistic distribution,
        pValue=mean(T<=T0);
    case 'right'
        % p=Pr(T>=t | H0) for a one-sided right-tail test-statistic distribution,
        pValue=mean(T>=T0);
    case 'both'
        % p=2*min{Pr(T>=t | H0),Pr(T>=t | H0)} for a two-sided test-statistic distribution. 
        % If the distribution of T is symmetric about zero, then p=Pr(abs(T)>=abs(t) | H0)

        %pValue=mean(abs(T-mean(T))>=abs(T0-mean(T)));
        pValue=2*min(mean(T<=T0),mean(T>=T0));
end

h=pValue<alpha;

%--------------------------------------------------------------------------

end

