function [ fSpeed ] = fSpeedSolve( Xt, Xv, model, mwpcaInitialize)
%FSPEEDSOLVE automatically determines the value of the eta parameter for
%RPCA and the window size parameter for MWPCA. This is done by first
%approximating the curve corresponding to the sum of the Q-statistics from
%a validation set for different parameter values, and then selecting the
%forgetting speed parameter value as that which obtains the largest
%orthogonal distance from the straight line drawn between the maximum and
%the minimum summed errors.
%
% Required input arguments: 
%   Xt    : A data matrix of observations in which rows represent observations, and columns represent variables.
%           Xt is used to obtain an initial PCA model using pcaModel.
%   Xv    : A data matrix used to validate parameter estimates.
%   model : An 'RPCA' or 'MWPCA' model structure to start monitoring from.
%
% Optional input arguments:
%   mwpcaInitialize : Function from mwpcaModel used to initialize an MWPCA
%                     model. 
%
% I/O: fSpeed = fSpeedSolve( Xt, Xv, model, []);    % For RPCA
% I/O: fSpeed = fSpeedSolve( Xt, Xv, model, @(Xt, fSpeed)mwpcaInitialize(Xt, fSpeed));    % For MWPCA
%
% The output of FSPEEDSOLVE is the value:
%   fSpeed :  The selected forgetting speed parameter eta or the window size for RPCA and MWPCA, respectively
%
% see also rpcaModel, mwpcaModel

[nv,p] = size(Xv);
midPoint = round(nv/2);
Xvt = Xv(1:midPoint, :);
Xvv = Xv((midPoint+1):end,:);
nXvv = size(Xvv,1);
[nt,p] = size(Xt);
if strcmp(model.type,'RPCA')
    upper = 0.9999;
%     upper = 1-2/(nt);
    lower = 0.9;
    tolX=1e-6;
    txt_label='2/(1-\eta)';
elseif strcmp(model.type,'MWPCA')
    %[nt,p] = size(Xt);
    upper = nt;
    lower = p+1;
    tolX=1;
    txt_label='H';
end

% minumum
[ fSpeed_min ] = fminbnd( @(x)fun_curve(x), lower, upper,optimset('TolX',tolX));

% % ortogonal distance
% Q1=[lower fun_curve(lower)];
% Q2=[upper fun_curve(upper)];
% [ fSpeed_ort ] = fminbnd( @(x)-fun_ort(x), lower, upper,optimset('TolX',tolX));

% grid
n_grid=30;
f=zeros(1,(n_grid+1));
if strcmp(model.type,'MWPCA')
    x=linspace(lower, upper, n_grid);
    for i=1:n_grid,
        f(i)=fun_curve(x(i));
    end
    f((n_grid+1))=fun_curve(fSpeed_min);
    fVec = [[x fSpeed_min]; f]';
    [~,I] = sort(fVec(:,1));
    fVec = fVec(I,:);
    fSpeed_min = fVec(find( fVec(:,2) == (min(fVec(:,2)))),1);
    fSpeed_minErr = (fVec(find( fVec(:,2) == (min(fVec(:,2)))),2));
    fSpeed_min = fSpeed_min(1);
    fSpeed_minErr = fSpeed_minErr(1);
    figure
    plot(fVec(:,1), fVec(:,2),'b',fSpeed_min, fSpeed_minErr,'ro')
    xlabel(txt_label);
    ylabel('Sum of Squared Preditiction Errors')
end

if strcmp(model.type,'RPCA')
    %x=linspace(log(lower), log(upper), n_grid);
    x=linspace(2/(1-(lower)), 2/(1-(upper)), n_grid);
    for i=1:n_grid,
        f(i)=fun_curve(1-2/(x(i)));
    end
    f((n_grid+1))=fun_curve(fSpeed_min);
    fVec = [[x 2/(1-fSpeed_min)]; f]';
    [~,I] = sort(fVec(:,1));
    fVec = fVec(I,:);
    fSpeed_min = 1-2/(fVec(find( fVec(:,2) == (min(fVec(:,2)))),1));
    fSpeed_minErr = (fVec(find( fVec(:,2) == (min(fVec(:,2)))),2));
    fSpeed_min = fSpeed_min(1);
    fSpeed_minErr = fSpeed_minErr(1);
    figure
    plot(fVec(:,1), fVec(:,2),'b',2/(1-fSpeed_min), fSpeed_minErr,'ro')
    xlabel(txt_label);
    ylabel('Sum of Squared Prediction Errors')

end

if strcmp(model.type,'MWPCA')
    fSpeed_min=round(fSpeed_min);
%     fSpeed_ort=round(fSpeed_ort);
end

if strcmp(model.type,'RPCA')
    fSpeed_min=fSpeed_min;
%     fSpeed_ort=fSpeed_ort;
end


% choice = menu('Parameter selection',['Minimum (',num2str(fSpeed_min),')'],['Ortogonal distance (',num2str(fSpeed_ort),')'],'Other (user defined)');
% switch choice
%     case 1
%         fSpeed=fSpeed_min;
%     case 2
%         fSpeed=fSpeed_ort;
%     case 3
%         name='Parameter selection';
%         prompt='Which value to use?';
%         numlines=1;
%         defaultanswer={''};
%         options.WindowStyle='normal';
%         answer=inputdlg(prompt,name,numlines,defaultanswer,options);
%         fSpeed=str2num(answer{1});
% end

choice = menu('Parameter selection',['Minimum (',num2str(fSpeed_min),')'],'Other (user defined)');
switch choice
    case 1
        fSpeed=fSpeed_min;
    case 2
        name='Parameter selection';
        prompt='Which value to use?';
        numlines=1;
        defaultanswer={''};
        options.WindowStyle='normal';
        answer=inputdlg(prompt,name,numlines,defaultanswer,options);
        fSpeed=str2num(answer{1});
end


    function [ f ] = fun_curve(fSpeed)
        if strcmp(model.type,'RPCA')
            model.eta = fSpeed;
            model = wPcaModel(Xt, model.eta, model.kSelect, model.kSelectType);
               
            model.aT = 0;
            model.aE = 0;
            model.ucl.sT = limitT( model.aT,model.k);
            model.ucl.sE = limitQ( model.L, model.k,  model.aE);
            
            [~, ~, model2] = rpcaMonitor( Xvt, model);
            [stat, ~] = rpcaMonitor(Xvv, model2);
        elseif strcmp(model.type,'MWPCA')
            [model] = mwpcaInitialize(Xt, fSpeed, model.cpvThresh);
            model.h = round(fSpeed);
            [~, ~, model2] = mwpcaMonitor( Xvt, model);
            [stat, ~] = mwpcaMonitor(Xvv, model2);
        end
        
        f = sum(stat.sE)/nXvv;
    end

    function [ f ] = fun_ort(fSpeed)
        P=[fSpeed fun_curve(fSpeed)];
        f= abs(det([Q2-Q1;P-Q1]))/norm(Q2-Q1); % for row vectors.
    end

end