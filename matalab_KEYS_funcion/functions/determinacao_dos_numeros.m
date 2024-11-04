function [N,HE] = determinacao_dos_numeros(id_valid,id_calibra, X , y)
%DETERMINACAO_DOS_NUMEROS Summary of this function goes here
%   Detailed explanation goes here

 %id_valid - matrix com as posiçoes dos dados que iram usar para validar (c , nº amostras para validar )
%  id_calibra - matrix com as posiçoes dos dados que iram usar 
% x - dados 
% y- resposta
tox_N=0.15;
C_max=length(id_valid(:,1));

for n=1:15
    for c=1:C_max
        C=c;
    idx=id_calibra(C,:);
    x_c=ones(length(idx),length(X(1,:)));
    Y_c=ones(length(idx),1);
    for i=1:length(idx)
        x_c(i,:)=X(idx(i),:);
        Y_c(i)=y(idx(i));
    end
    idx=id_valid(C,:);
    x_v=ones(length(idx),length(X(1,:)));
    y_v=ones(length(idx),1);
    for i=1:length(idx)
        x_v(i,:)=X(idx(i),:);
        y_v(i)=y(idx(i));
    end
    N=n;
    [ ~,  ~, ~, ~ , BETA, ~ ] = plsModel( x_c, Y_c, N );
    
    XS_V=x_v * BETA.T; 
    XS_C=x_c* BETA.T;
    model_classMDL = fitcdiscr(XS_C,Y_c,'CrossVal', 'off');
%     model_classMDL = fitcdiscr(XS_C,Y_c,'CrossVal', 'off','discrimType','pseudoLinear');
    [predictedLabels]=predict(model_classMDL, XS_V); 
    HE(c,n)=mean(y_v == predictedLabels);


% F1 score 
    % Create confusion matrix
%     cm = confusionmat(y_v , predictedLabels);
%     cmt=cm';
%     dig=diag(cmt);
%     sun_rows= sum(cmt,2);
%     sum_col=sum(cmt,1);
% 
%     precision= mean(dig ./ sun_rows); %determiar a precião
%     recall=mean(dig ./ sun_rows);
% %     
%     HE(c,n) = 2 * ((precision*recall)/(precision + recall));

    
    end
  
    if N>1
        p=signrank(HE(:,N-1),HE(:,N),'tail','right');
        if p<tox_N
            N=N-1;
            break
        end
    end



% F1 score 
% % Create confusion matrix
% cm = confusionmat(trueLabels, predictedLabels);
% cmt=cm';
% dig=diag(cmt);
% sun_rows= sum(cmt,2);
% sum_col=sum(cmt,1);
% 
% precision= mean(dig ./ sun_rows); %determiar a precião
% recall=mean(dig ./ sun_rows);
% 
% F1(c,n) = 2 * ((precision*recall)/(precision + recall));
end

