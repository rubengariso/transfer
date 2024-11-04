function [X] = multiblock(x,method)

switch method
    case 'NS'
       KB=1;
    case 'SBS'
        n=size(x,2);
        KB=1/(n^(0.25));
    case 'HBS'
        n=size(x,2);
        KB=1/(n^(0.5));
    case 'SHBS'
        n=size(x,2);
        KB=1/(n);
    case 'SNVS'
        n=size(x,2);
        sigma=var(x);
        den=sum(sigma.^2);
        KB=(n^(0.25))/den^0.5;
    case 'HBVS'
        sigma=var(x);
        den=sum(sigma.^2);
        KB=1/den^0.5;
     case 'SHBVS'
        sigma=var(x);
        den=sum(sigma.^2);
        KB=1/den;
    otherwise
        fprintf('Wrong data pretreat method!')
      
end
X=x/KB;
end