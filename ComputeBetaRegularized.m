function [ Beta ] = ComputeBetaRegularized(X,Y,eps)
% Find the optimal parameters beta in the multilinear
% regression Y = B'*X + Epsilon
% (I. Horenko, 2015)

% Compute the covariance matrices sum(Xt*Xt') and sum(Xt*Yt')
[n,k]=size(X);
XX=zeros(k,k);YY=zeros(size(Y,2),size(Y,2));
for t=1:n
    XX=XX+X(t,:)'*X(t,:);
    YY=YY+X(t,:)'*Y(t,:);

end
% Fit the Multilinear Regression
Beta=inv(XX+eps^2*eye(k))*YY;

end

