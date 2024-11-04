function [ stat, ucl, model ] = maxDMonitor( X, model )
%MAXDMONITOR Summary of this function goes here
%   Detailed explanation goes here

% Yeh et al. (2005) -------------------------------------------------------

S0z=model.S0z;% start as S0=eye(p);
lambda=model.lambda;
AnZ=model.AnZ;% iS=inv(Sigma); diS=diag(diag(iS)); AnZ=diS^-0.5*iS; 

[n,m]=size(X);

stat.maxD=zeros(n,1);
ucl.maxD=model.ucl.maxD*ones(n,1);
for i=1:n,
    
    z=AnZ*X(i,:)';
    Sz=lambda*(z*z')+(1-lambda)*S0z;
    
    D1=sum((diag(Sz)-1).^2);
    D2=sum((vmat(Sz)-0).^2);
    
    stat.maxD(i)=max((D1-model.mewD1)/model.sigmaD1,(D2-model.mewD2)/model.sigmaD2);
    
    % update values
    S0z=Sz;
end

% update model
model.S0z=Sz;


%--------------------------------------------------------------------------

    function [ r ] = vmat( R )
                
        nx=size(R,1);
        C=ones(nx,nx);
        C=tril(C,-1);
        [k, j]=find(C==1);
        n0=[j k];
        nmax=size(n0,1);
        % ordem 0
        r=zeros(nmax,1);
        for j=1:nmax,
            r(j)=R(n0(j,1),n0(j,2));
        end
        
end

end

