function [F,RMSE,Yhat] = piecewise_direct_standardization_v2(X_new,X_ref,windoh,lambda_ref,lambda_new)


% 29-10-24

%% make sure that the \lambda start in the lower and end in the lagger value

% [lambda_ref,idx]=sort(lambda_ref);
% X_ref=X_ref(:,idx); 
% [lambda_new,idx]=sort(lambda_new);
% X_new=X_new(:,idx);

n_ref=length(lambda_ref);
n_new=length(lambda_new);
resolution_new=abs(mean(diff(lambda_new)));
window=ceil(windoh/resolution_new);

scaling_method='mean-centering';
F=zeros(n_new,n_ref);
% Faux=cell(n_ref,1);
for i=1:n_ref
    % lambda_begin_interval=lambda_ref(i)-windoh;
    % lambda_end_interval=lambda_ref(i)+windoh;
    % [~,j_min] = min(abs(lambda_new-lambda_begin_interval));
    % [~,j_max] = min(abs(lambda_new-lambda_end_interval));

    [~,index] = min(abs(lambda_new-lambda_ref(i)));
    index_lb= max(index-window,1);
    index_ub=min(index+window,n_new);

    % try
    [ kpls ] = kSelectPLS_MCCV_3d(X_new(:,index_lb:index_ub),X_ref(:,i),[],scaling_method );
    % catch
    %     gg=o;
    % end
    % meanY=mean(X_ref);
    % stdY=ones(1,size(X_ref,2));
    % zY=zScale(X_ref,meanY,stdY);
    [ ~,~,~,~,BETA, ~ ] = plsModel(X_new(:,index_lb:index_ub), X_ref(:,i) ,kpls );
    F(index_lb:index_ub,i)=BETA.Y;
    % Faux{i}=BETA.Y;
end

% for i=1:n_ref,
%     [~,index] = min(abs(lambda_new-lambda_ref(i)));
%     index_lb= max(index-window,1);
%     index_ub=min(index+window,n_new);
%     F(index_lb:index_ub,i)=Faux{i};
% end

Yhat=X_new*F;
E=X_ref-Yhat;
RMSE=sqrt(mean(E.^2,2));

end

