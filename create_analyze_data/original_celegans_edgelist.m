clear all
close all
clc

load celegans277.mat
W = celegans277matrix;


i=0
iter = 1;
attempts = 1;
myIter = i;
R = W;
eff = 0
fname = ['edge_list_randomized_myiter_', num2str(myIter), '_attempts_', num2str(attempts), '_effRewirings_',num2str(eff), '.csv']
fid = fopen(fname, 'w');
[a,b] = find(R>0);
for r = 1: size(a,1)
      fprintf(fid, '%d %d\n', a(r), b(r));
end
fclose(fid)
