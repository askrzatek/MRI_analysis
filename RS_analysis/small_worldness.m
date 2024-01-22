function sw = small_worldness(matric_bin,rep)


% Calculate Path length  divide pathlength of random network (same
% node-edge

% matric_bin : le graphe binarisé



% get random network with same size of nodes and edges

K = length(find(matric_bin == 1))/2;
N = length(matric_bin);

C = clustering_coef_bu(matric_bin);
D = distance_bin(matric_bin);
infinite_dist=false
[lambda,efficiency,ecc,radius,diameter] = charpath(D,0,infinite_dist);

for i =1:50
randNetwork = makerandCIJ_und(N,K);
Cr = clustering_coef_bu(randNetwork);
Dr = distance_bin(randNetwork);
[lambdar] = charpath(Dr,0,infinite_dist);
sw_rep(i)= (mean(C)/mean(Cr))/(lambda/lambdar)

end
sw_ = sw_rep(~isinf(sw_rep)); 
sw = mean(sw_)

end