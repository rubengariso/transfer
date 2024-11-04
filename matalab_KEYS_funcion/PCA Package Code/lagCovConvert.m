
function [ S, vec, lag, p, m, covmat ] = lagCovConvert( covmat, Lmax, L )

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

end

in_cov=[true(1,m) false(1,p-m)];
for i=m+1:p,
    if lag(i)<=L(vec(i)),
        in_cov(i)=true;
    end
end

covmat=covmat(in_cov,in_cov);
S=eig(covmat);
S=sort(S,'descend');

end
