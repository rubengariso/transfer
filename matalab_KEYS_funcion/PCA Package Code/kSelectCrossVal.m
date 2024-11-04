function [ k ] = kSelectCrossVal( X, mode )
% KSELECTCROSSVAL selects the number of latent variables to be used on PCA 
% models base on cross-validation.
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
%   X  : A data matrix of observations where rows represent observations, and columns represent variables. 
% mode : Method used to determine the PRESS, wheiter 5-fold ('Kfold') or apriximation('aprox'). 
%
% I/O: [ k ] = kSelectCrossVal( X )
%
% The output of KSELECTCROSSVAL is:
%   k : The selected number of latent variables

[n,m]=size(X);

if nargin<=1,
    mode='aprox';
end

if strcmp(mode,'aprox')==1,
    vecL = svd(cov(X));
    PRESS0=sum(vecL);
    PRESS=zeros(m,1);
    for i=1:m,
        PRESS(i)=sum(vecL(i+1:end));% proportional to PRESS
    end
    
elseif strcmp(mode,'Kfold')==1,
    
    [ PRESS, SE ] = pressCrossVal( X, 5 );
    PRESS=PRESS(:);
    PRESS0=SE(1);
    
end

NFI=zeros(m,1);
P=zeros(m,1);
for i=1:m,
    NFI(i)=1-PRESS(i)/PRESS0;
    DR=(n-i-1)*(m-i);
    Q=DR/(n*m-m);
    P(i)=NFI(i)*Q;
end

[~,k]=max(P);

    function [ PRESS, SE ] = pressCrossVal( X, k )
        
        [iter, nv]=size(X);
        kit=floor(iter/k);
        
        EX=X(1:kit*k,:);
        R=zeros(1,nv);
        PRESS=zeros(1,nv);
        SE=zeros(1,nv);
        for ii=1:nv,
            
            PRESSpar=zeros(1,k);
            for j=1:k,
                
                EXtrain=zeros(kit*(k-1),nv);
                g=1;
                for r=1:k,
                    if r~=j,
                        EXtrain((g-1)*kit+1:g*kit,:)=EX((r-1)*kit+1:r*kit,:);
                        g=g+1;
                    else
                        EXtest=EX((r-1)*kit+1:r*kit,:);
                    end
                end
                
                [ V, ~ ]=eig(cov(EXtrain));
                V=V(:,end);
                E=EXtest*(eye(nv)-V*V');
                PRESSpar(j)=sum(sum(E.^2));
            end
            PRESS(ii)=sum(PRESSpar);
            SE(ii)=sum(sum(EX.^2));
            R(ii)=PRESS(ii)/SE(ii);
            
            
            [ V, ~ ]=eig(cov(EX));
            V=V(:,end);
            EX=EX*(eye(nv)-V*V');
            
        end
        
        M=size(EX,1);
        PRESS=PRESS/M;
        SE=SE/M;
                
    end

end

