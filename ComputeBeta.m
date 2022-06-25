function [ Beta ] = ComputeBeta(X,Y)
% Find the optimal parameters beta in the multilinear
% regression Y = B'*X + Epsilon
% (I. Horenko, 2015)

% Compute the covariance matrices sum(Xt*Xt') and sum(Xt*Yt')
XX=zeros(size(X,2),size(X,2));YY=zeros(size(Y,2),size(Y,2));
n=size(X,1);
Y
for t=1:n
    XX=XX+X(t,:)'*X(t,:);
    YY=YY+X(t,:)'*Y(t,:);
    %YY=cov(X(t:1)',Y(t:1));
end
%XX
disp("AND NOW YYYYY")
YY
% Fit the Multilinear Regression
Beta=inv(XX)*YY;

end

