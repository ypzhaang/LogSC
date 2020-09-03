clear all;
clc;
rng('default'); rng(1);
%load('myMnist_50.mat');
load('pltYaleB_50.mat');

param.lambda1=0.1;
param.lambda2=0.01;
param.lambda3=1;
param.lambda4=1;
param.nAtom=128;

profile on
[ D,X ] = LogSC(trData,param);
profile viewer

save('Yale_X.mat','X');