fea = rand(50,70);
gnd = [ones(10,1);ones(15,1)*2;ones(10,1)*3;ones(15,1)*4];

res = kmeans(fea,4);

res = bestMap(gnd,res);
%=============  evaluate AC: accuracy ==============
AC = length(find(gnd == res))/length(gnd);
%=============  evaluate MIhat: nomalized mutual information =================
MIhat = MutualInfo(gnd,res);