function output = TDseries(n,m,p,T_params,A_params,scenario,fault)
%%%% Description
    % This function simulates time-dependent series of data. 5 settings are
    % possible. 
        %'static': A static series with no time dynamics
        %'AR': A series with autocorrelation
        %'ARI': A series with an ARI(1,1) process
        %'NSS': A series with a non-stationary model for the loadings
        %'ARI_NSS': A series combining ARI and NSS.
        
%%%% INPUTS
%%%%%%%% Inputs for scenario
%%%%AR
    %% T_params.phi: vector of AR coefficients
%%%%ARI
    %% T_params.phi: see above
    %% T_params.mu: vector of mu values for the ARI model
%%%%NSS
    %% A_params.freq: vector of frequencies to use as input for the rotation matrix
%%%%ARI_NSS
    %% T_params.phi: see above
    %% T_params.mu: vector of mu values for the ARI model
    %% A_params.freq: vector of frequencies to use as input for the rotation matrix
%%%%%%%% Inputs for faults
    % A fault can be introduced to one of the scores. It's type and
    % specifications are defined by the user in the 'fault' input. 
    % fault.scenario: the scenario of the fault:
        %'step': increase in the mean of one of the scores
        %'ramp': slow increase over time
        %'scale': change the importance of one of the scores so that the
        %         correlation relationship changes
        %'dynamic': change on the dynamic of one of the scores
        %'correlation': adds correlation between two scores
        %'error': adds error with a diferent signal-to-noise ratio to one
        %         of the sensors
    % fault.loc: The observation corresponding to the location where the fault begins.
    % fault.shift: The magnitude of the final shift in the step and ramp in
    % standard deviations
    % fault.speed: how many observations it takes for the ramp fault to
    % fault.rescale: how much the scale of the first score should change if
    % the scenario is 'scale'.
    % fault.phi: new AR coefficient for one of the scores.
    % fault.rho: new correlation between two scores 
    % fault.SNR: new signal-to-noise ratio (default SNR=20 dB)
    % reach the final fault.shift value.
    % fault.type: which level of the process the fault takes place on
        %'score': fault occurs at the score level
        %'sensor': fault occurs on a variable in the final data
    % fault.refSigma: standard deviation to normalize the fault's magnitude (default: 'error'):
        %'signal': total variance of the signal
        %'error': variance of the error
        

% Change log
% 28.05.2013 - Changes on T0 are scaled by their standard desviation.
%            - sX was replaced by the theoritical value for 'static', 'AR' and 'ARI'
% 02.07.2013 - Added a 'dynamic' change to one of the scores.
% 04.07.2013 - Added a 'correlation' change on the scores.
%            - Added a 'error' change on X
% 18.07.2013 - The reference sT and sX are stored in T_params.base
% 04.03.2015 - Faults are now scaled based on the error standard deviation
% 17.03.2015 - Corrected the scaling procedure introduced previously.
%            - Unnecessary commands are now a comments.
    
if nargin == 6
    fault.loc = n;
    fault.type = '';
    fault.refSigma = 'error';
end
%This part is for the scores

if isfield(T_params,'base')==0,
    sT=[];
    sX=[];
else
    sT=T_params.base.stdT;
    sX=T_params.base.stdX;
end
sigma_eps=0.1;
    
    switch scenario
        case {'static', 'NSS'} 
            T0=sigma_eps*random('normal',0,1,n,p);
            if isempty(sT)==1,
                sT=sigma_eps*ones(1,p);
            end
        case {'AR','AR_NSS'} 
            phi=T_params.phi;
            % generate scores
            T0=zeros(n,p);
            Et=0.1*random('normal',0,1,n,p);
            for i=2:n,
                T0(i,:)=phi.*T0(i-1,:)+Et(i,:);
            end
            if isempty(sT)==1,
                sT=sqrt(sigma_eps^2./(1-phi.^2));
            end
         case {'ARI','ARI_NSS'}
            phi=T_params.phi;
            mu=T_params.mu;
            
            % generate scores
            T0=zeros(n,p);
            Et=0.1*random('normal',0,1,n,p);
            for i=3:n,
                T0(i,:)=mu+T0(i-1,:)+phi.*(T0(i-1,:)-T0(i-2,:))+Et(i,:);
            end
            if isempty(sT)==1,
                sT=std(T0);
            end
        case {'IMA'}
            phi=T_params.phi;
            mu=T_params.mu;
            
            % generate scores
            T0=zeros(n,p);
            Et=0.1*random('normal',0,1,n,p);
            for i=2:n,
                T0(i,:)=mu+T0(i-1,:)+Et(i,:)-phi.*Et(i-1,:);
            end
            if isempty(sT)==1,
                sT=std(T0);
            end
        case {'MA'}
            phi=T_params.phi;
            
            % generate scores
            T0=zeros(n,p);
            Et=0.1*random('normal',0,1,n,p);
            for i=2:n,
                T0(i,:)=Et(i,:)-phi.*Et(i-1,:);
            end
            if isempty(sT)==1,
                sT=sqrt((1+phi.^2)*sigma_eps^2);
            end
    end
    %save base standard deviation of T
    T_params.base.stdT=sT;
    
    %This part adds faults to the first score
    if strcmp(fault.refSigma,'signal')==1,
        sigmaRefT=sT(1);
    else
        sigmaRefT=sigma_eps;
    end
    if strcmp(fault.type,'score')
        switch fault.scenario
            case 'step'
                T0(fault.loc:end,1)=T0(fault.loc:end,1)+fault.shift*sigmaRefT;
            case 'ramp'
                increase = linspace(0,fault.shift*sigmaRefT, fault.speed);
                T0(fault.loc:(fault.loc+fault.speed-1),1) = T0(fault.loc:(fault.loc+fault.speed-1),1)+increase';
                T0((fault.loc+fault.speed):end,1) = T0((fault.loc+fault.speed):end,1)+fault.shift*sigmaRefT;
            case 'scale'
                T0(fault.loc:end,1) = T0(fault.loc:end,1)*fault.rescale;
            case 'dynamic'
                for i=fault.loc:n,
                    T0(i,1)=fault.phi.*T0(i-1,1)+Et(i,1);
                end
            case 'correlation'
                sigma=sqrt(1+fault.rho.^2);
                R=[1 fault.rho]'./sigma;
                T0(fault.loc:end,1)=T0(fault.loc:end,1:2)*R;
        end
    end
    
%This part is for the A matrix
    
    if isfield(A_params,'base')==0,
        A=random('normal',0,1,m,m);
        [U,S,V]=svd(A);
        P0=V(:,1:p);
        A_params.base=P0;
    else
        P0=A_params.base;
    end

    switch scenario 
        case {'static','AR','ARI','IMA','MA'}
             Xr=T0*P0';
             if isempty(sX)==1 && strcmp(scenario,'ARI')==0 && strcmp(scenario,'IMA')==0,
                 S=P0*diag(sT.^2)*P0';
                 sX=sqrt(diag(S))';
             elseif isempty(sX)==1,
%                  SNR=20*log10(1/0.1);
%                  Xnoise=awgn(Xr,SNR,'measured');
%                  sX=std(Xnoise-Xr)/0.1;% 0.1 is a scaling factor
                 sX=std(Xr);
%                     S=P0*diag(sT.^2)*P0';
%                     sX=sqrt(diag(S))';
             end
        case {'NSS','AR_NSS','ARI_NSS'}    
            % rotate model
            A_params.freq = 0.001;
            n_ang=length(A_params.freq);
            t=0:n-1;
            theta=(15/180)*pi*sin(2*pi*repmat(t',1,n_ang).*repmat(A_params.freq,n,1));
            Xr=zeros(n,m);
             for i=1:n,
                R=eye(m);
                for j=1:n_ang,
                    r=[cos(theta(i,j)) -sin(theta(i,j)); sin(theta(i,j)) cos(theta(i,j))];
                    Rk=eye(m);
                    Rk(j:j+1,j:j+1)=r;
                    R=R*Rk;
                end
                Pk=R*P0;
                Xr(i,:)=T0(i,:)*Pk';
             end
             if isempty(sX)==1,
                 sX=std(Xr(1:fault.loc,:));
             end
    end
    
    % Add noise to X
    frac_E=0.23;% E has 5% of total variance in iid
%     frac_E=0.001;% E has 5% of total variance in iid
    sigma_eps_X=sqrt(diag(P0*sigma_eps.^2*eye(p)*P0'));
    sigma_E=frac_E*mean(sigma_eps_X);
%     sigma_total_X=sqrt(sT.^2+sigma_E^2);
    
    
    E=sigma_E*random('normal',0,1,n,m);
    Xr=Xr+E;
    %save base standard deviation of X
    T_params.base.stdX=sX;
    
    
    if strcmp(fault.refSigma,'signal')==1,
        sigmaRefX=sX(1);
    else
        sigmaRefX=sigma_E;
    end
    if strcmp(fault.type,'sensor')
        switch fault.scenario
            case 'step'
                Xr(fault.loc:end,1)=Xr(fault.loc:end,1)+fault.shift*sigmaRefX;
            case 'ramp'
                increase = linspace(0,fault.shift*sigmaRefX, fault.speed);
                Xr(fault.loc:(fault.loc+fault.speed-1),1) = Xr(fault.loc:(fault.loc+fault.speed-1),1)+increase';
                Xr((fault.loc+fault.speed):end,1) = Xr((fault.loc+fault.speed):end,1)+fault.shift*sigmaRefX;
%             case 'error'
%                 sE=sX(1)/10^(fault.SNR/20);
%                 Xr(fault.loc:end,1) = Xr(fault.loc:end,1)-E(fault.loc:end,1)+sE*random('normal',0,1,n-fault.loc+1,1);% removes the default error and adds new one
        end
    end
    
    
    output.X = Xr;
    output.T = T0;
    output.T_params = T_params;
    output.A_params = A_params;