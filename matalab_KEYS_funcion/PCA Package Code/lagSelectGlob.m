function [ L, KSV, KSVR ] = lagSelectGlob( X, Lmax, mode )
%LAGSELECT uses a method of sucessive SVD decompositins to select the maximum
% number of lags more approriate for the data. 
%
% Reference: 
%
%    Rato, T. J. and M. S. Reis (2013),
%    "Defining the structure of DPCA models and its impact on process monitoring and prediction activities.",
%    Chemometrics and Intelligent Laboratory Systems, 125, 74-86.
%
% Required input arguments:
%   X         : Data matrix (observations X variables)
%   Lmax      : Maximum number of lags.
%   nit       : Maximum number of iterations. nit also defines the total of lagged
%               variables added. To perform the full selection, set nit as empty.
%   mode      : Mode to select the number of lags:
%               1: Automatic;
%               2: Manual.
%
% I/O: [ L, LSV, LSVR ] = lagSelectGlob( X, Lmax, mode );
%
% The output of LAGSELECTGLOB are:
%   L    : Selected number of lags.
%   KSV  : Key singular value.
%   KSVR : Key singular value ratio.
%
% Example:
% To determine the number of lags to use on a data matrix X, use:
% Lmax=20;% up to 20 lags.
% mode='manual';% plot the optimization function.
% [ L ] = lagSelectGlob( X, Lmax, nfig )
%
% see also lagSelect.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin==1 || isempty(Lmax)==1,
    Lmax=10;
end

if nargin<3,
    mode='auto';
end

m=size(X,2);

KSV=zeros(Lmax+1,1);
KSVR=[ NaN; zeros(Lmax,1)];
c=NaN*ones(Lmax+1,1);

[ Z ] = lagData(X, Lmax );
covmat=cov(Z);
[ ~, info ] = lagCov( covmat, Lmax*ones(1,m), Lmax*ones(1,m) );
for j=0:Lmax,
    
    [ S ] = lagCov( [], [], j*ones(1,m), info );
    KSV(j+1)=S(end-m+1);
    if j>=1,
        KSVR(j+1)=KSV(j+1)/KSV(j);
    end
    
    if j>=2 && KSVR(j+1)<KSVR(j),
        c(j+1)=1;
    end
    
end

dto1=(KSV-min(KSV))./(max(KSV)-min(KSV));
dto2=(KSVR-min(KSVR))./(max(KSVR)-min(KSVR));
w=[1 1];
DTO=sqrt(w(1)*dto1.^2+w(2)*dto2.^2);
DTO2=DTO.*c;
[~, L]=min(DTO2);
L=L-1;% because it starts on 0

% graphics

if strcmp(mode,'manual')==1,
    
    figure
    subplot(1,2,1);
    [AX,H1,H2]=plotyy(0:Lmax,KSV,0:Lmax,KSVR);
    set(H1,'LineStyle','--','Marker','.')
    set(H2,'LineStyle',':','Marker','o')
    xlabel('lag')
    set(get(AX(1),'Ylabel'),'String','KSV')
    set(get(AX(2),'Ylabel'),'String','KSVR')
%     set(AX,{'ycolor'},{'k';'k'})
    subplot(1,2,2);
    plot(1:Lmax,DTO(2:end),'.-',1:Lmax,DTO2(2:end),'ro')
    xlabel('Stage')
    ylabel('\phi')
    
    prompt={'Which stage to use?'};
    name='Select number of lags';
    numlines=1;
    defaultanswer={''};
    options.WindowStyle='normal';
    answer=inputdlg(prompt,name,numlines,defaultanswer,options);
    L=str2double(answer{1});
    
end

function [ S, info ] = lagCov( covmat, Lmax, L, info )

if nargin<=3,
    d=length(L);
    p=size(covmat,1);
    vec=lagData(repmat(1:d,max(Lmax)+1,1),Lmax);
    lag=zeros(1,p);
    ct=d;
    for i=1:max(Lmax),
        lag(ct+1:ct+sum(Lmax>=i))=i;
        ct=ct+sum(Lmax>=i);
    end
    info.vec=vec;
    info.lag=lag;
    info.p=p;
    info.m=d;
    info.covmat=covmat;
end

in_cov=[true(1,info.m) false(1,info.p-info.m)];
for i=info.m+1:info.p,
    if info.lag(i)<=L(info.vec(i)),
        in_cov(i)=true;
    end
end

covmat=info.covmat(in_cov,in_cov);
S=eig(covmat);
S=sort(S,'descend');

end

end



