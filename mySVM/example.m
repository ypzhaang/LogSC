clear all;
clc;

load iris;
libsvmwrite('svm_iris.dat',dataLabs,sparse(dataIns));

clear all;
[tag,fea]=libsvmread('svm_iris.dat');

SVMModel=libsvmtrain(tag,fea,'-s 0 -t 2');
predict=libsvmpredict(tag,fea,SVMModel);

error=sum(predict~=tag)/length(tag)