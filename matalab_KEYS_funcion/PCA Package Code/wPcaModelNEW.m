function [model] = wPcaModelNEW(X, eta, kSelect, kSelectType, varargin)
    [n,m] = size(X);     
    
    %initialize using the first observation
    %The initially centered observation is at exactly zero
    xi_center = zeros(1,m);
    %The initial observation's values are the average
    model.loc = X(1,:);
    %There is no variation
    model.S = zeros(m);
    model.eta = eta;
    model.kSelect = kSelect;
    model.kSelectType = kSelectType;
    model.n = 1;
    model.cpvThresh = 0.95;
    if nargin>1,
        varargin=reshape(varargin,2,[]);
        for c=1:size(varargin,2),
            switch varargin{1,c}
                case 'cpvThresh'
                    model.cpvThresh=varargin{2,c};
            end
        end
    end

%This loop builds up the correct weighted covariance matrix.
    for i = 2:n,
            eta_n = model.n/(model.n+1)*model.eta; % weight parameter
            xi = X(i,:); % in LV, this will be the signal from the new seal
            xi_center = xi - model.loc;% mean centering of the new signal

            model.loc_old = model.loc; % store the old mean vector b
            model.loc = eta_n*model.loc + (1-eta_n)*xi; % update b (this should be done using the raw data, not a centered version)
            d_loc = model.loc - model.loc_old; % difference between new and old b
            model.S = eta_n*(model.S+d_loc'*d_loc)+(1-eta_n)*xi_center'*xi_center; % update for the covariance matrix to be applied for scaling
            %model.R = corrcov(model.S);
            model.n = model.n + 1; % update of the number of samples
 
    end
                model.R = corrcov(model.S);
                model_i = pcaModelNEW(model.R, model.kSelect, model.kSelectType, model.cpvThresh);
                model.P = model_i.P;
                model.L = model_i.L;
                model.k = model_i.k;
                model.type = 'RPCA'; 
                 

end