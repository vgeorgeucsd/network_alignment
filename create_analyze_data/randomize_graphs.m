clear all
close all
clc

load celegans277.mat
W = celegans277matrix;


for i=[0,1,10,100,1000,10000,100000]
iter = 1;
attempts = 1;
myIter = i;
[R, eff] = randmio_dir_connected(W,iter,myIter,attempts);
eff
fname = ['edge_list_randomized_myiter_', num2str(myIter), '_attempts_', num2str(attempts), '_effRewirings_',num2str(eff), '.csv']
fid = fopen(fname, 'w');
[a,b] = find(R>0);
for r = 1: size(a,1)
      fprintf(fid, '%d %d\n', a(r), b(r));
end
fclose(fid)
end
