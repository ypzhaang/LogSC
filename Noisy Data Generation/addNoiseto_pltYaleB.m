
clear all
clc;

load pltYaleB_50;

XXXX=[];
for i=1:1:500
    X=reshape(trData(:,i),[32,32]);
    %imagesc(X); colormap(gray);
    X=X/255;
    XX=imnoise(X,'salt & pepper',0.1);%gaussian
    XX=XX*255;
    %imagesc(XX); colormap(gray);
    XXX=reshape(XX,[1,32*32]);
    XXXX=[XXXX;XXX];
end
clear trData;
trData=XXXX';
save('pltYaleB_50_sp.mat','trData');