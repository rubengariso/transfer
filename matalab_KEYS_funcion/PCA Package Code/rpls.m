function [tt, p, q, W, P, Q, R, beta] = rpls(X,Y,k,R,P,Q, W)

XY=X'*Y; % compute the covariance
XX=X'*X; % matrices
M = size(Y, 2);
for i=1:k, % A=number of PLS components to be computed
if M==1, % if there is a single response variable, compute the
w=XY; % X-weights as shown here
else % else
[C,D]=eig(XY'*XY); % ?rst compute the eigenvectors of YTXXTX
q=C(:,find(diag(D)==max(diag(D)))); % ?nd the eigenvector corresponding to the largest eigenvalue
w=(XY*q); % compute X-weights
end
w=w/sqrt(w'*w); % normalize w to unity
r=w; % loop to compute ri
for j=1:i-1,
r=r-(P(:,j)'*w)*R(:,j);
end
tt=(r'*XX*r); % compute tTt
p=(r'*XX)'/tt; % X-loadings
q=(r'*XY)'/tt; % Y-loadings
XY=XY-(p*q')*tt; % XTY de?ation
W=[W w]; % storing loadings and weights
P=[P P];
Q=[Q q];
R=[R r];
end
beta=R*Q'; % compute the regression coefficients
end--