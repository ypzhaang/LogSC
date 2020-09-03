%function [ output_args ] = LogSC(Y,Lambda)
%LRGSR Summary of this function goes here
%   Detailed explanation goes here
% This function is to solve the following problem:
%||Y-DX||_2;F+lambda1*sum||x_i||_1
%+0.5*lambda2*sum(W_ij*||x_i-x_j||_2;2)
%+lambda3*||Z||_*+lambda4*||E||_1
%subject to:
%            W=Z'+Z; Y=YZ+E;  Z>=0; d_i=1 for all i;
%input:
%         Y: the sample matrix m*n, where m is the dimensional and n is
%         the number of samples;
%         D: the dictionary for sparse representation
%         X: the sparse representation matrix of Y on D
%output:
%         D: the learned dictionary.
% Written by Yupei Zhang on 2017/05
% Contact : xjtuiezhyp@stu.xjtu.edu.cn
function [ D,X ] = LogSC(trData,param)

Y=normc(trData);
%% settings
lambda1=param.lambda1;
lambda2=param.lambda2;
lambda3=param.lambda3;
lambda4=param.lambda4;
nAtom=param.nAtom;

[M,N]=size(Y); % M-dimensional N samples

tol1=1e-4;   % threshold for the error in constraint
tol2 = 1e-5; %threshold for the change in the solutions
max_mu = 1e6;
mu = min(M,N)*tol2;
step_mu=1.9;
norm2Y = norm(Y,2);
normfY=norm(Y,'fro');
%eta needs to be larger than ||X||_2^2, but need not be too large.
eta = norm2Y*norm2Y*1.02;

maxIter=1000;
DIS=1;
%% Initializings
D=rand(M,nAtom); % random dictionary
D=normc(D);
LassoParm.lambda=lambda1;
X=mexLasso(Y,D,LassoParm);
E = sparse(M,N);
Z=zeros(N,N);%eye(N,N);%
W=Z'+Z;

M1=zeros(N,N);
M2=zeros(M,N);

%% Learning
t=0;
while t<maxIter
    t=t+1;
    Zt=Z; Et=E; Dt=D; Xt=X; Wt=W;
    %1 update E******************************************************
    QE=Y-Y*Z+M2/mu;
    %E=max(0,QE-lambda4/mu)+min(0,QE+lambda4/mu);%E_1
    E=solve_l1l2(QE,lambda4/mu);%E_l2,1
    %2 update Z*******************************************************
    QZ_1=Y-E+M2/mu;
    QZ_2=W+M1/mu;
    Grad=2*Z+2*Z'-QZ_2-QZ_2'+Y'*Y*Z-Y'*QZ_1;

    ZG=Z-Grad/eta;
    [U,S,V]=svd(ZG,'econ');
    S=diag(S);
    svp=length(find(S>lambda3/(mu*eta)));
    if svp>=1
        S=S(1:svp)-1/(mu*eta);
    else
        svp=1;
        S=0;
    end
    Z=U(:,1:svp)*diag(S)*V(:,1:svp)';
    %Z=max(0,Z);
    %3 update W*******************************************************
    QW=Z'+Z-M1/mu;
    %B=distMat(X',X');
    B=squareform(pdist(X'));
    W=(2*mu*QW-lambda2*B)/(2*mu);
    W=max(0,W);
    %imagesc(W);colormap(flipud(gray));
    %4 update X********************************************************
        % unnormalized Laplacian;
        DCol = full(sum(W,2));
        DL = spdiags(DCol,0,speye(size(W,1)));
        L = DL - W; 
        %...
        XX=[];
        for i=1:1:N
            y=Y(:,i);
            Qx=X*L(i,:)'-L(i,i)*X(:,i);
            y_new=[y;sqrt(lambda2*L(i,i))*Qx];
            D_new=[D;-sqrt(lambda2*L(i,i))*eye(size(Qx,1),nAtom)];
            x=mexLasso(y_new,D_new,LassoParm);
            XX=[XX,x];
        end
        X=XX;
    %5 update D******************************************************
    D=learn_basis(Y,X,1);
    
    %<<<<<<<<checking the stoping>>>>>>>>>>>>>>>>>>>>
    relChgZ=norm(Zt-Z,'fro')/normfY;
    relChgE=norm(Et-E,'fro')/normfY;
    relChgZE=max(relChgZ,relChgE);
    
    relChgD=norm(Dt-D,'fro')/normfY;
    relChgX=norm(Xt-X,'fro')/normfY;
    relChgDX=max(relChgD,relChgX);
    
    %relChgW=norm(Wt-W,'fro')/normfY;
    relChg=max(relChgDX,relChgZE);
    %relChg=max(max(relChgDX,relChgZE),relChgW);
    
    
    errLRR=norm(Y-Y*Z-E,'fro')/normfY;
    errDict=norm(Y-D*X,'fro')/normfY;
    err=max(errLRR,errDict);
    
    convergenced= relChg<tol1 && err<tol1;
    
    if DIS
        if t==1 || mod(t,10)==0 || convergenced
        disp(['t=' num2str(t) ',E=' num2str(norm(E,'fro')) ...
            ',rankZ=' num2str(svp) ...
            ',relChg=' num2str(relChg) ',err=' num2str(err)]);
        end
    end
    %<<<checking end>>>
    
    if convergenced
        break;
    else   
        %6 update M************************************************************
        M1=M1+mu*(W-Z'-Z);
        M2=M2+mu*(Y-Y*Z-E);
        %7 update parameters**************************************************
        if mu*relChgZE<tol2
            mu=min(max_mu,mu*step_mu);
        end
    end    
end

end
