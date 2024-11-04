function [model] = iplsBackwardModel_class(X, Y, rndIndex, nInter)
% 14.04.2023

threshold=0.10;

numVar=size(X,2);

nRCV=length(rndIndex);

breakPoints=round(linspace(1,numVar+1,nInter+1));

varIndex=1:numVar;
intervalIndex=nan(1,numVar);
for i=1:nInter
    intervalIndex(varIndex>=breakPoints(i) & varIndex<breakPoints(i+1))=i;
end

% Initial model -----------------------------------------------------------

[ model_tent, HR_tent, HR_weigted_tent ] = plsModel_class(X ,Y, rndIndex);

model_ref=model_tent;
HR_ref=HR_tent;
HR_weigted_ref=HR_weigted_tent;
med_HR=median(HR_weigted_ref);    

%--------------------------------------------------------------------------

% Backward stepwise -------------------------------------------------------

inModel=true(1,numVar);
removedModel_intervals=[];
for i=1:nInter-1,
    
    pVal=nan(nInter,1);
    HR_tent=nan(nRCV,nInter);
    HR_weigted_tent=nan(nRCV,nInter);
    medHR=nan(nInter,1);
    model_tent=cell(nInter,1);
    
    k_range=setdiff(1:nInter,removedModel_intervals);
    parfor k=1:nInter,
        if sum(k_range==k)==1,
            inModel_tent=inModel;
            inModel_tent(intervalIndex==k)=false;
            
            [ model_tent{k}, HR_tent(:,k), HR_weigted_tent(:,k) ] = plsModel_class(X(:,inModel_tent) ,Y, rndIndex);
            
%          [P,H] = signrank(...,'tail',TAIL) performs the test against the
%          alternative hypothesis specified by TAIL:
%           'both'  -- "median is not zero (or M)" (two-tailed test, default)
%           'right' -- "median is greater than zero (or M)" (right-tailed test)
%           'left'  -- "median is less than zero (or M)" (left-tailed test)
            % left: H1: HR_weigted_ref < HR_weigted_tent(:,k)
            pVal(k) = signrank(HR_weigted_ref,HR_weigted_tent(:,k),'tail','left');% For a two-sample test, the alternate hypothesis states the data in x – y come from a distribution with median greater than zero.
            medHR(k)=median(HR_weigted_tent(:,k));
        end
    end
    
    % interval can be removed from the model if HR_weigted_ref < HR_weigted_tent(:,k):
    % that is, H0 is rejected:
    can_remove=pVal<threshold;
    if sum(can_remove)>=1,
        % remove the interval that leades to the best improvement
        medHR(~can_remove)=-inf;
        [~,indexSelect]=max(medHR);
        
        model_ref=model_tent{indexSelect};
        HR_ref=HR_tent(:,indexSelect);
        HR_weigted_ref=HR_weigted_tent(:,indexSelect);
        
        inModel(intervalIndex==indexSelect)=false;
        removedModel_intervals=[removedModel_intervals indexSelect];
        
    else
        break
    end
    
    if i==nInter-1,
       warning('attempted to remove all interval!') 
    end
end

model=model_ref;
model.inModel=inModel;

%--------------------------------------------------------------------------

end

