function [model, HR] = plsModel_Regress(X, Y, rndIndex)
% 28.09.2023
%ruben 
kmax=size(X,2);

nRCV=length(rndIndex);

thresh_kpls=0.15;

hr=nan(nRCV,kmax);
hr_weigted=nan(nRCV,kmax);
for k=1:kmax,
    for r=1:nRCV,
        stdY=std(Y(rndIndex{r}.indexTrain));
        meanY=mean(Y(rndIndex{r}.indexTrain));
        zy=(Y(rndIndex{r}.indexTrain)-meanY)/stdY;
        [ ~, ~, ~,~, BETA ] = plsModel( X(rndIndex{r}.indexTrain,:), zy, k );
        Yhat =X(rndIndex{r}.indexTest,:)*BETA.Y*stdY+meanY;
   
        % 
        % 
        % 
        % 
        % [ P, Q, B, W, BETA ] = plsModel( X(rndIndex{r}.indexTrain,:), Y(rndIndex{r}.indexTrain), k );
        % 
        % T=X*BETA.T;
        % 
        % classMDL = fitcdiscr(T(rndIndex{r}.indexTrain,:),Y(rndIndex{r}.indexTrain),'CrossVal','off');
        % 
        % Yhat = predict(classMDL,T(rndIndex{r}.indexTest,:));
      %% Root Mean Squared Error
        hr(r,k) =sqrt( mean(Y(rndIndex{r}.indexTest)-Yhat).^2);

        % TP=Y(rndIndex{r}.indexTest)==Yhat & Y(rndIndex{r}.indexTest)==1;
        % TN=Y(rndIndex{r}.indexTest)==Yhat & Y(rndIndex{r}.indexTest)==0;
        % hr_weigted(r,k) = (1/2)*(sum(TP)/sum(Y(rndIndex{r}.indexTest)==1)+sum(TN)/(sum(Y(rndIndex{r}.indexTest)==0)));

    end
    
    if k>1,
        
        %          [P,H] = signrank(...,'tail',TAIL) performs the test against the
%          alternative hypothesis specified by TAIL:
%           'both'  -- "median is not zero (or M)" (two-tailed test, default)
%           'right' -- "median is greater than zero (or M)" (right-tailed test)
%           'left'  -- "median is less than zero (or M)" (left-tailed test)
            % left: H1: hr_weigted(:,k-1) < hr_weigted(:,k)
            % right: H1: hr_weigted(:,k-1) > hr_weigted(:,k)
            
% %         pVal= signrank(hr(:,k),hr(:,k-1),'tail','left');% For a two-sample test, the alternate hypothesis states the data in x – y come from a distribution with median less than 0.

%         % right:
%         % H0: hr_weigted(:,k-1) <= hr_weigted(:,k)
%         % H1: hr_weigted(:,k-1) > hr_weigted(:,k)
%         pVal= signrank(hr_weigted(:,k-1),hr_weigted(:,k),'tail','right');% For a two-sample test, the alternate hypothesis states the data in x – y come from a distribution with median greater than zero.
%         % if H0 is rejected, stop ; rejected means => hr_weigted(:,k-1) > hr_weigted(:,k)
%         % if H0 is not rejected, continue for k+1 ; not rejected means => hr_weigted(:,k-1) <= hr_weigted(:,k)
%         if pVal<thresh_kpls,
%             kpls=k-1;
%             break
%         end
        pVal= signrank(hr_weigted(:,k-1),hr_weigted(:,k),'tail','right');
        % left:
        % H0: hr_weigted(:,k-1) >= hr_weigted(:,k)
        % H1: hr_weigted(:,k-1) < hr_weigted(:,k)
        % pVal= signrank(hr_weigted(:,k-1),hr_weigted(:,k),'tail','left');
        % if H0 is rejected, continue for k+1 ; rejected means => hr_weigted(:,k-1) < hr_weigted(:,k)
        % if H0 is not rejected, stop ; not rejected means => hr_weigted(:,k-1) >= hr_weigted(:,k)
        if pVal>thresh_kpls,
            kpls=k-1;
            break
        end
    end
    
end

% train final model
stdY=std(Y);
meanY=mean(Y);
zy=(Y(rndIndex{r}.indexTrain)-meanY)/stdY;
[ P, Q, B, W, BETA ] = plsModel( X, zy, kpls );

% [ P, Q, B, W, BETA ] = plsModel( X, Y, kpls );
% T=X*BETA.T;
% classMDL = fitcdiscr(T,Y,'CrossVal','off');

model.k=kpls;
model.P=P;
model.Q=Q;
model.B=B;
model.W=W;
model.BETA=BETA;
% model.classMDL=classMDL;
model.inModel=true(1,kmax);

HR=hr(:,kpls);
% HR_weigted=hr_weigted(:,kpls);

end

