clear all;
clc;

load pltYaleB_50.mat;

Y=trData;

imagesc(Y'*Y); colormap(jet);
set(gca,'yticklabel',[]);
set(gca,'xticklabel',[]);
box on;
%colorbar;

