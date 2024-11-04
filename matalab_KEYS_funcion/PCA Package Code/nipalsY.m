function [W,P,Q,R,beta]=nipalsY(X,Y,A) % A=number of PLS latent vectors to be computed
% beta=PLS regression coef?cients
convcrit=1*exp(1)^-10; % convergence criteria for the NIPALS algorithm
for i=1:A,
u=Y(:,find(std(Y)==max(std(Y)))); % assigning the y-variables with the largest variance
% to u_initial
diff=1;
while (diff>convcrit),
u0 = u;
w=X'*u/(u'*u); % compute w
w=w/(sqrt(sum(w.^2))); % normalizing w to unity
r=w; % compute r via equation (30) so that the score
if i == 1,
    W = []; Q = []; R = []; T = []; U = []; P = [];
end
    
if i>1, % vector t can be computed directly from the
for j=1:i-1, % original X
r = r-(P(:,j)'*w)*R(:,j);
end
end
t=X*r; % compute t
q=Y'*t/(t'*t); % compute q
u=Y*q/(q'*q); % compute u
diff=(sqrt(sum((u0-u).^2)))/(sqrt(sum(u.^2)));
end
p=X'*t/(t'*t); % compute X-loadings
W=[W w]; % storing loadings and weights
R=[R r];
P=[P p];
Q=[Q q];
T=[T t];
U=[U u];
Y=Y-t*q'; % only Y block is de?ated
end
beta=R*Q'; % compute the regression coefficients
end