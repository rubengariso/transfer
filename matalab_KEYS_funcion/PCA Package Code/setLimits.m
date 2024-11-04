function [ UCLp, alpha_p, alfaout ] = setLimits( stat, aG )
% SETLIMITS selects the control limits for the set of monitoring statistics
% in stat (observations x statistics) in order to get the desired false
% detection rate (aG), subject to an equal false detection rate for each 
% monitoring statistic.
%
% Required input arguments:
%      stat : Monitoring statistics (observations x statistic)
%        aG : Target false detection rate
%
% I/O: [ UCLp, alfap, alfaout ] = setLimits( stat, aG );
%
% The output of SETLIMITS are:
%      UCLp : UCL limits for each monitoring statistic (1 x statistic) 
%    alfa_p : False detection rate of each monitoring statistic (1 x
%             statistic)
%   alfaout : Global false detection rate obtained (1 x statistic)



[nit,nst]=size(stat);
alpha=fzero(@(x)F(x),aG/nst);

[ Fm, UCLp, alpha_p ] = F ( alpha );
alfaout=Fm+aG;

    function [ f, UCLp, alpha_p ] = F ( alpha )
        
        if alpha<eps,
            alpha=eps;
        end
        
        UCLp=zeros(1,nst);
        alpha_p=zeros(1,nst);
        v=zeros(nit,nst);
        for s=1:nst,
            UCLp(s)=prctile(stat(:,s),(1-alpha)*100);
            v(:,s)=stat(:,s)>UCLp(s);
            alpha_p(s)=sum(v(:,s))/nit;
        end
        
        v=max(v,[],2);
        t=sum(v)/nit;
        f=t-aG;
        
    end


end

