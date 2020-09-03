clear all;
clc;

load COIL20.mat;

trData=[];
for i=1:1:10
    D=fea(1+(i-1)*72:i*72,:);
    trData=[trData; D(1:50,:)];
end
trData=trData';
save('pltCOIL20.mat','trData');