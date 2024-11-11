function [F,RMSE,Yhat] = piecewise_direct_standardization(X_source,X_destination,lambda_source,lambda_destination,windoh,index_mccv)


% 29-10-24

%% make sure that the \lambda start in the lower and end in the lagger value

% [lambda_ref,idx]=sort(lambda_ref);
% X_ref=X_ref(:,idx); 
% [lambda_new,idx]=sort(lambda_new);
% X_new=X_new(:,idx);

n_destination=length(lambda_destination);
n_source=length(lambda_source);
resolution_new=abs(mean(diff(lambda_source)));
window=ceil(windoh/resolution_new);

scaling_method='mean-centering';
F=zeros(n_source,n_destination);
% Faux=cell(n_ref,1);
SE=zeros(size(index_mccv.train,2),n_destination);
n_SE=zeros(size(index_mccv.train,2),n_destination);
for i=1:n_destination
    % lambda_begin_interval=lambda_ref(i)-windoh;
    % lambda_end_interval=lambda_ref(i)+windoh;
    % [~,j_min] = min(abs(lambda_new-lambda_begin_interval));
    % [~,j_max] = min(abs(lambda_new-lambda_end_interval));

    [~,index] = min(abs(lambda_source-lambda_destination(i)));
    index_lb= max(index-window,1);
    index_ub=min(index+window,n_source);

    % try
    [ kpls, ~, SE_k, n_SE_k  ] = kSelectPLS_MCCV_3d(X_source(:,index_lb:index_ub),X_destination(:,i),index_mccv,scaling_method );
    SE(:,i)=SE_k(:,end);
    n_SE(:,i)=n_SE_k(:,end);
    % catch
    %     gg=o;
    % end
    % meanY=mean(X_ref);
    % stdY=ones(1,size(X_ref,2));
    % zY=zScale(X_ref,meanY,stdY);
    [ ~,~,~,~,BETA, ~ ] = plsModel(X_source(:,index_lb:index_ub), X_destination(:,i) ,kpls );
    F(index_lb:index_ub,i)=BETA.Y;
    % Faux{i}=BETA.Y;
end

% for i=1:n_ref,
%     [~,index] = min(abs(lambda_new-lambda_ref(i)));
%     index_lb= max(index-window,1);
%     index_ub=min(index+window,n_new);
%     F(index_lb:index_ub,i)=Faux{i};
% end

Yhat=X_source*F;
% E=X_ref-Yhat;
% RMSE=sqrt(mean(E.^2,2));
RMSE=sqrt(sum(SE,2)./sum(n_SE,2)); 
end

