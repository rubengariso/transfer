function [ NL, fmin ] = lagSelect( X, m_outputs, Lmax, nit, mode )

%LAGSELECT uses a method of sucessive SVD decompositins to select the number 
% of lags more approriate for each variable in X. A different number of lags
% is selected for each variable
%
% Reference: 
%
%    Rato, T. J. and M. S. Reis (2013),
%    "Defining the structure of DPCA models and its impact on process monitoring and prediction activities.",
%    Chemometrics and Intelligent Laboratory Systems, 125, 74-86.
%
% Required input arguments:
%   X         : Data matrix (observations X variables)
%   m_outputs : Number of inputs (on the first columns of X, i.e. X=[Xoutputs Xinputs]).
%               If the number of outputs is unknown define m_outputs as empty.
%   Lmax      : Maximum number of lags for all variables.
%   nit       : Maximum number of iterations. nit also defines the total of lagged.
%               variables added. To perform the full selection, set nit as empty.
%   mode      : Mode to select the number of lags:
%               1: Automatic;
%               2: Manual.
%
% I/O: [ NL, fmin ] = lagSelect( X, m_outputs, Lmax, nit, mode );
%
% The output of LAGSELECT are:
%   NL   : Vector containing the number of lags for each variables.
%   fmin : Smallest singular value of a SVD decomposition with the lags
%          defines by NL.
%
% Example:
% To determine the number of lags to use on a data matrix X, use:
% m_outputs=[];% no distinction between inputs and outputs.
% Lmax=10;% up to 10 lags for each variables.
% nit=[];% the algorithm only stops when Lmax lags are added to all variables
% mode='manual';% the functions plots the key singular values, key singular values
% ratios and optimization function. The user is prompt to select the most
% suitable stage (usually the minimum of the optimization function).
% [ NL, fmin ] = lagSelect( X, m_outputs, Lmax, nit, mode )
%
% see also lagSelectGlob.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin==4,
    mode='auto';
end
w=[1 1];

if length(Lmax)==1,
    Lmax=Lmax*[1 1];
end

if isempty(m_outputs)==1,
    m_outputs=size(X,2);
end

if isempty(nit)==1,
    nit=max(Lmax)*size(X,2)+10;
end

LMx=Lmax(1);
LMu=Lmax(2);

imax=size(X,3);
NLx=zeros(imax,m_outputs);
f=zeros(imax,1);
for i=1:imax,
    [ NLx(i,:), f(i) ] = searchLag( X(:,1:m_outputs,i), LMx, nit, mode, 0, 0, w );
end

[~,ind]=min(f);
NLx=NLx(ind,:);

if size(X,2)>m_outputs,
    
    Lmaxu=[NLx repmat(LMu,1,size(X,2)-m_outputs)];
    
    NL=zeros(imax,size(X,2));
    f=zeros(imax,1);
    for i=1:imax,
        [ NL(i,:), f(i) ] = searchLag( X(:,:,i), Lmaxu, nit, mode, NLx, 1, w );
    end
    [fmin,ind]=min(f);
    NL=NL(ind,:);
    
else
    
    NL=NLx;
    fmin=f(ind);
    
end

end

function [ L0, f0 ] = searchLag( Xref, Lmax, imax, mode, nx, isu, w )

if length(nx)>1 || isu==1,
    lag=nx;
    nx=length(nx);
else
    lag=0;
end

% pre-process data ------------------------------------------------

% auto-scaling
m=size(Xref,2);
[ Z ] = zscore( Xref );

%------------------------------------------------------------------

% determine minimum singular value of SVD -------------------------

if isu==0,
    Lm=zeros(1,m);
    Lm(1:nx)=lag;
    Lm(nx+1:end)=Lmax;
    Lmax=Lm;
    
else
    Lm=Lmax;
end
[ ZLmax ] = lagData( Z, Lm );
covmat=cov(ZLmax);
[ S, info ] = lagCov( covmat, Lmax, Lm );
Smin=S(end);

%------------------------------------------------------------------

% determine individual lags ---------------------------------------

L0=zeros(1,m);
L0(1:nx)=lag;

h=zeros(1,m);
h(1:nx)=inf;
s=zeros(1,imax);
r=zeros(1,imax);
c=NaN*ones(1,imax);
Ls=zeros(imax,m);
caux=0;

for i=1:imax,
    
    % KSV and KSVR of correct stage
    if i==1,
        [ S ] = lagCov( [], [], L0, info );
        s(i)=S(end);
    else
        s(i)=Sind(ind);
    end
    Ls(i,:)=L0;
    
    if i==1,
        r(i)=NaN;
    elseif i==2,
        r(i)=s(i)/s(i-1);
    else
        r(i)=s(i)/s(i-1);
        if r(i)<r(i-1) || caux==1;
            c(i)=1;
            caux=1;
        end
    end
    
    if abs(s(i)-Smin)<=0,
        break
    end
    
    % singular value for each added lag
    Sind=zeros(1,m);
    for j=1:m,
        
        if h(j)~=inf,
            K=zeros(1,m);
            K(j)=1;
            L1=L0+K;
            
            [ S ] = lagCov( [], [], L1, info );
            
            Sind(j)=S(end);
            
            h(j)=Sind(j);
            
        end
        
    end
    
    % select variable to add lag
    [~,ind]=min(h);
    L0(ind)=L0(ind)+1;
    
    % avoid setting more lags than the maximum
    if L0(ind)==Lmax(ind),
        h(ind)=inf;
    end
    
end
imax=i;

%------------------------------------------------------------------

opt1=(s-min(s))./(max(s)-min(s));
opt2=(r-min(r))./(max(r)-min(r));
w=w/sum(w);
OPT1=sqrt(w(1)*opt1.^2+w(2)*opt2.^2);
OPT2=sqrt(w(1)*opt1.^2+w(2)*opt2.^2).*c;


if strcmp(mode,'auto')==1,
    [~,ind]=min(OPT2);
    L0=Ls(ind,:);
    f0=s(ind);
end

% manual selection ------------------------------------------------

if strcmp(mode,'manual')==1,
    s=s(1:imax);
    r=r(1:imax);
    figure
    subplot(1,2,1)
    [AX,H1,H2]=plotyy(0:imax-1,s,0:imax-1,r);
    set(H1,'LineStyle','--','Marker','.')
    set(H2,'LineStyle',':','Marker','o')
    xlabel('Stage')
    set(get(AX(1),'Ylabel'),'String','KSV')
    set(get(AX(2),'Ylabel'),'String','KSVR')
    % set(AX,{'ycolor'},{'k';'k'})% axis collor
    legend('KSV','KSVR')
    
    subplot(1,2,2)
    plot(0:imax-1,OPT1(1:imax),':.b',0:imax-1,OPT2(1:imax),'.r');
    xlabel('Stage')
    ylabel('\phi');
    
%     ind=input('Which stage to use?','s');
%     ind=str2double(ind)+1;
    
    prompt={'Which stage to use?'};
    name='Select number of lags';
    numlines=1;
    defaultanswer={''};
    options.WindowStyle='normal';
    answer=inputdlg(prompt,name,numlines,defaultanswer,options);
    ind=str2double(answer{1})+1;
    
    L0=Ls(ind(:),:);
    f0=s(ind(:));
    
end

%------------------------------------------------------------------

end

function [ S, info ] = lagCov( covmat, Lmax, L, info )

if nargin<=3,
    m=length(L);
    p=size(covmat,1);
    vec=lagData(repmat(1:m,max(Lmax)+1,1),Lmax);
    lag=zeros(1,p);
    c=m;
    for i=1:max(Lmax),
        lag(c+1:c+sum(Lmax>=i))=i;
        c=c+sum(Lmax>=i);
    end
    info.vec=vec;
    info.lag=lag;
    info.p=p;
    info.m=m;
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
