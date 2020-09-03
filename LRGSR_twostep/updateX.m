function X=updateX(Y,D,X_k,L,lambda1,lambda2)
%UPDATEX Summary of this function goes here
%   Detailed explanation goes here
% this function is to solve the following problem:
%min ||Y-DX||_2;F+lambda1*sum||x_i||_1+lambda2*Tr(XLX')


[m,n]=size(Y);
X=[];
I=eye(size(X_k,1),size(D,2));
param=[];
param.mode=2;
param.lambda=lambda1;
param.lambda2=0;
for i=1:1:n
    Q=X_k*L(i,:)'-X_k(:,i)*L(i,i);
    Q=Q/L(i,i); 
    if isempty(Q) || L(i,i)==0
        Y_i_e=Y(:,i);
        D_e=D;
    else
        Y_i_e=[Y(:,i);sqrt(lambda2*L(i,i))*Q];
        D_e=[D;-sqrt(lambda2*L(i,i))*I];
    end
    
    Y_i_e=myNormalization(Y_i_e,'column');
    D_e=myNormalization(D_e,'column');
    X_i=mexLasso(Y_i_e,D_e,param); 
    %X_i=mexOMP(Y_i_e,D_e,param);
    
    X=[X,X_i];
end

end

