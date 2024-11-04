function [model] = iplsForwardModel_class_mod_V2(X, Y, rndIndex,nInter,intervalIndex)
% 19.01.2024

threshold=0.10;

kmax=size(X,2);

nRCV=length(rndIndex);

% breakPoints=round(linspace(1,kmax+1,nInter+1));
% 
% varIndex=1:kmax;
% intervalIndex=nan(1,kmax);
% for i=1:nInter
%     intervalIndex(varIndex>=breakPoints(i) & varIndex<breakPoints(i+1))=i;
% end

% Initial model -----------------------------------------------------------

for i=1:nInter,
    
   inModel_tent=intervalIndex==i;
   
    [ model_tent, HR_tent, HR_weigted_tent ] = plsModel_class_mod(X(:,inModel_tent) ,Y, rndIndex);
    
    if i==1,
        model_ref=model_tent;
        HR_ref=HR_tent;
        HR_weigted_ref=HR_weigted_tent;
        med_HR=median(HR_weigted_ref);
        int_ref=1;
    else
        if med_HR<median(HR_weigted_tent)
            model_ref=model_tent;
            HR_ref=HR_tent;
            HR_weigted_ref=HR_weigted_tent;
            med_HR=median(HR_weigted_ref);
            int_ref=i;
        end
    end
    
end

%--------------------------------------------------------------------------

% Forward stepwise --------------------------------------------------------

inModel=intervalIndex==int_ref;
inModel_var=int_ref;
for i=1:nInter,
    
    pVal=nan(nInter,1);
    HR_tent=nan(nRCV,nInter);
    HR_weigted_tent=nan(nRCV,nInter);
    medHR=nan(nInter,1);
    model_tent=cell(nInter,1);
    k_range=setdiff(1:nInter,inModel_var);
    parfor k=1:nInter,
        if sum(k_range==k)==1,
            inModel_tent=inModel | intervalIndex==k;
            
            [ model_tent{k}, HR_tent(:,k), HR_weigted_tent(:,k) ] = plsModel_class_mod(X(:,inModel_tent) ,Y, rndIndex);
            
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
    
    % interval can enter in the model if HR_weigted_ref < HR_weigted_tent(:,k):
    % that is, H0 is rejected:
    can_enter=pVal<threshold;
    if sum(can_enter)>=1,
        % add the best
        medHR(~can_enter)=-inf;
        [~,indexSelect]=max(medHR);
        
        model_ref=model_tent{indexSelect};
        HR_ref=HR_tent(:,indexSelect);
        HR_weigted_ref=HR_weigted_tent(:,indexSelect);
        
        inModel=inModel | intervalIndex==indexSelect;
        inModel_var=[inModel_var indexSelect];
        
    else
        break
    end
    
end

model=model_ref;
model.inModel=inModel;

%--------------------------------------------------------------------------

end

