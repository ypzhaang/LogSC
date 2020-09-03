clear all;
clc;

load('pltYaleB_50.mat');
%load('myMnist_50.mat');

Y=normc(trData);

lambda=2;
[Z,E] = ladmp_lrr(Y,lambda);
W=max(Z+Z',0);
%imagesc(W);

nAtom=256;
alpha=0.1;
beta=0.01;
nIter=100;
Binit=Y(:,1:nAtom);
[D X stat] = GraphSC(Y, W, nAtom, alpha, beta, nIter, Binit);
save('YaleB_LogSC2.mat','D','X');