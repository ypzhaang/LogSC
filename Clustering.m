%clear all;
%clc;

path(path,'./LogSC');
path(path,'./ClusterEV');
path(path,'./Data');
%load YaleB_X.mat;%pltYaleB_50;%

%load Yale.mat;
load COIL20.mat

%fea = X';
gnd=gnd;
%gnd=[];
%for i=1:1:10
%    gnd=[gnd;ones(50,1)*i];
%end

AC=[];MI=[];
for t=1:1:10
res = kmeans(fea,2,'Distance','cosine');%,'Distance','cosine'
res = bestMap(gnd,res);
ACt = length(find(gnd == res))/length(gnd);
AC=[AC, ACt];
MIhat = MutualInfo(gnd,res);
MI=[MI, MIhat];
end

avgAC=mean(AC);
avgMI=mean(MI);
disp(['avgAC=' num2str(avgAC) ', avgMIhat=' num2str(avgMI)]);