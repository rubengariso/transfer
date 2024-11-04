clear
clc
close all

DO_model=0;

%% data --------------------------------------------------------------------

n_test_noc=500;
n_test_noc=500;
n_test_fault=500;

index_noc=[true(n_test_noc,1); false(n_test_fault,1) ];
index_fault=[false(n_test_noc,1); true(n_test_fault,1) ];
scenario='ARI';
    
% file=['C:\Users\EricS\Documents\MATLAB\Research\ProcessSimulation\TDseries\',scenario,'\',scenario,'(R).mat'];
file=['C:\Users\mebios\Documents\TDseries\',scenario,'\',scenario,'(R).mat'];
load(file,'data_ref','data_val')

Xref=[data_ref.X];
Xval=data_val.X;
aG = 0.01;
[mwpcaModelInit] = mwpcaModel(Xref(1:500,:,1), Xref(501:2000,:,1), aG, 'h', 396);
[rpcaModelInit] = rpcaModel(Xref(1:500,:), Xref(501:2000,:), aG, 'lambda', 0.995);


%==========================================================================

%% Faults =================================================================

model={	'PCA'            'T^2_P_C_A'                 'Q_P_C_A';...
        'DPCA_LS1'       'T^2_D_P_C_A_,_L_S_1'       'Q_D_P_C_A_,_L_S_1';...
        'DPCA_LS2'       'T^2_D_P_C_A_,_L_S_2'       'Q_D_P_C_A_,_L_S_2';...
        'DPCADR'         'T^2_p_r_e_v'               'T^2_r_e_c';...
        'RPCA'           'T^2_R_P_C_A'               'Q_R_P_C_A';...
        'MWPCA'           'T^2_M_W_P_C_A'            'Q_M_W_P_C_A'};

n_model=size(model,1);
DO_plot=0;    
nfig=1;


ID_faults=[0:6 10:14];
for ID=ID_faults,

    
%     file=['C:\Users\EricS\Documents\MATLAB\Research\ProcessSimulation\TDseries\',scenario,'\',scenario,'(',num2str(ID),').mat'];
    file=['C:\Users\mebios\Documents\TDseries\',scenario,'\',scenario,'(',num2str(ID),').mat'];
    load(file,'data_fault')
    load(file,'data_start')
    
    % pre-processing ----------------------------------------------------------
    
    X=data_fault.X;
    Xref = data_start.X;
   [nit,~,rep]=size(X);
   % RPCA model -----------------------------------------------------------
    UCL.RPCA = NaN(nit,2,rep);
    if isempty(rpcaModelInit)==0,
         STAT.RPCA=nan(nit,2,rep);
        for r=1:rep,
            [RPCA_model] = rpcaModel(Xref(500:1000,:,r), Xref(1001:2000,:,r), aG, 'lambda', rpcaModelInit.lambda, 'aT', rpcaModelInit.aT,'aE', rpcaModelInit.aE );
            [~, ~, RPCA_model] = rpcaMonitor(Xref(2001:end,:,r), RPCA_model); 
            [stat, ucl, ~] = rpcaMonitor(X(:,:,r), RPCA_model);
            STAT.RPCA(:,1,r) = stat.sT;
            STAT.RPCA(:,2,r) = stat.sE;
            UCL.RPCA(:,1,r) = ucl.uclT;
            UCL.RPCA(:,2,r) = ucl.uclE;
        end
    end
% 
%     %--------------------------------------------------------------------------
%     % MWPCA model -----------------------------------------------------------

    UCL.MWPCA = NaN(nit,2,rep);
    if isempty(mwpcaModelInit)==0,
         STAT.MWPCA=nan(nit,2,rep);
        for r=1:rep,
        [MWPCA_model] = mwpcaModel(Xref(200:1000,:,r), Xref(1001:2000,:,r), aG, 'h', mwpcaModelInit.h, 'aT', mwpcaModelInit.aT,'aE', mwpcaModelInit.aE);
        [stat, ucl, MWPCA_model] = mwpcaMonitor(Xref(2001:end,:,r), MWPCA_model); 
        [stat, ucl,~] = mwpcaMonitor(X(:,:,r), MWPCA_model);
        STAT.MWPCA(:,1,r) = stat.sT;
        STAT.MWPCA(:,2,r) = stat.sE;
        UCL.MWPCA(:,1,r) = ucl.uclT;
        UCL.MWPCA(:,2,r) = ucl.uclE;
        end
    end
  
    %--------------------------------------------------------------------------

    
    % TDR and FDR
    for i=5:n_model,
        TDR.(model{i,1})=nan(rep,2);
        FDR.(model{i,1})=nan(rep,2);
    end
    
    for r=1:rep,
        for i=5:n_model,
            if (i==3 && isempty(DPCA_LS1_model)==1) || (i==4 && isempty(DPCADR_model)==1),
                break
            end
            S1=STAT.(model{i,1})(:,1,r)./UCL.(model{i,1})(:,1,r);
            S2=STAT.(model{i,1})(:,2,r)./UCL.(model{i,1})(:,2,r);
            
            FDR.(model{i,1})(r,1)=sum(S1(index_noc)>=1)/sum(index_noc);
            TDR.(model{i,1})(r,1)=sum(S1(index_fault)>=1)/sum(index_fault);
            FDR.(model{i,1})(r,2)=sum(S2(index_noc)>=1)/sum(index_noc);
            TDR.(model{i,1})(r,2)=sum(S2(index_fault)>=1)/sum(index_fault);
            
            if DO_plot==1 && r==1,
                figure(nfig)
                nfig=nfig+1;
                subplot(1,2,1)
                plot(1:length(S1),S1,1:length(S2),S2)
                hold on
                plot(xlim,[1 1],':k')
                plot(n_test_noc*[1 1],ylim,':k')
                hold off
                title(['TDR(',model{i,2},')=',num2str(TDR.(model{i,1})(r,1)),', TDR(',model{i,3},')=',num2str(TDR.(model{i,1})(r,2))])
                legend(model{i,2},model{i,3})
                subplot(2,2,2)
                autocorr(S1(isnan(S1)==0 & index_noc))
                title(model{i,2})
                subplot(2,2,4)
                autocorr(S2(isnan(S2)==0 & index_noc))
                title(model{i,3})
            end
        end
    end
    
    % ROC
    for i=5:n_model,
        AUC.(model{i,1})=nan(rep,2);
    end
    for r=1:rep;
        
        for i=5:n_model,
            if (i==3 && isempty(DPCA_LS1_model)==1) || (i==4 && isempty(DPCADR_model)==1),
                break
            end
            
            if DO_plot==1 && r==1,
                figure(nfig)
                nfig=nfig+1;
                subplot(1,2,1)
                [x_roc,y_roc,~,AUC.(model{i,1})(r,1)] = perfcurve(index_fault,STAT.(model{i,1})(:,1,r),true(1));
                plot(x_roc,y_roc)
                hold on
                plot([0 1], [0 1],':k')
                hold off
                xlabel('False positive rate');
                ylabel('True positive rate')
                title(model{i,2})
                subplot(1,2,2)
                [x_roc,y_roc,~,AUC.(model{i,1})(r,2)] = perfcurve(index_fault,STAT.(model{i,1})(:,2,r),true(1));
                plot(x_roc,y_roc)
                hold on
                plot([0 1], [0 1],':k')
                hold off
                xlabel('False positive rate');
                ylabel('True positive rate')
                title(model{i,3})
            else
                [~,~,~,AUC.(model{i,1})(r,1)] = perfcurve(index_fault,STAT.(model{i,1})(:,1,r),true(1));
                [~,~,~,AUC.(model{i,1})(r,2)] = perfcurve(index_fault,STAT.(model{i,1})(:,2,r),true(1));
            end
        end
    end
    
    % save results
    
    file=['C:\Users\EricS\Documents\MATLAB\Research\ProcessSimulation\TDseries\',scenario,'\',scenario,'(',num2str(ID),')_STAT_EricAuto.mat'];
    save(file,'STAT','model','UCL','TDR','FDR','AUC')
    
end