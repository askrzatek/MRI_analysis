function zx = mean_connection_strength(z)


nodes = length(z);
zx= (sum(sum(z)))/((nodes*(nodes-1))/2);


end