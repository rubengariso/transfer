function [ T, model ] = missingData( X, model, MDmethod )
%MISSINGDATA applys missing data estimation to the first m variables in X
% based on a PCA model.
%
% Reference: 
%
%   Nelson, P. R. C., P. A. Taylor, et al. (1996),
%   "Missing data methods in PCA and PLS: Score calculations with incomplete observations.",
%   Chemometrics and Intelligent Laboratory Systems 35, 45-65.
%
% Required input arguments: 
%   X     : A data matrix of observations where rows represent observations, and columns represent variables. 
%   model : A structure witht the DPCA model parameters.
%
% I/O: [ T, model ] = missingData( X, model );
%
% The output of MISSINGDATA are:
%   T     : Estimated scores.
%   model : Orignal input model with a model.MDmatrix field correspondent
%           to the missing data model. Stored to save computation.

if nargin<=2,
   MDmethod='CMR'; 
end

if isfield(model,'MDmatrix')==0,
    
    switch MDmethod
        case 'CMR'
            %CMR
            Pp=model.P(model.m+1:end,:);
            a=size(model.L,1);% K - total number of variables
            B=[eye(model.k) zeros(model.k,a-model.k)];
            MDmatrix=B*model.L*Pp'*inv(Pp*model.L*Pp');
            
        case 'PMP'
            %PMP
            Pp=model.P(model.m+1:end,1:model.k);
            MDmatrix=inv(Pp'*Pp)*Pp';
            
        case 'TSR'
            %TSR
            P=model.P(model.m+1:end,:);
            Pp=model.P(model.m+1:end,1:model.k);
            L=model.L;
            Lk=model.L(1:model.k,1:model.k);
            MDmatrix=Lk*Pp'*Pp*inv(Pp'*P*L*P'*Pp)*Pp';
            
        case 'PLS'
            
            if model.k~=0
                Y=X*model.P(:,1:model.k);
                Z=X(:,model.m+1:end);
                [ ~, ~, ~, ~, BETA ] = fitpls( Z, Y, 10 );
                MDmatrix=BETA.Y';
            else
                MDmatrix=zeros([],size(model.L,1)-model.m);
            end
            
        case 'PLS-X'
            
            Y=X(:,1:model.m);
            Z=X(:,model.m+1:end);
            
            % select number of LV
            KF=5;
            [ R_PLS ] = CrossValPLS( Z, Y, KF );
            [ ~, ~, ~, ~, BETA ] = fitpls( Z, Y, R_PLS );
            MDmatrix=([BETA.Y eye(size(Z,2))]*model.P(:,1:model.k))';
            model.BETA=BETA;
            
    end
    
    model.MDmatrix=MDmatrix;
      
end

if isempty(X)==0,
    Zp=X(:,model.m+1:end);
    
    T=model.MDmatrix*Zp';
    T=T';
else
    T=[];
end

end


