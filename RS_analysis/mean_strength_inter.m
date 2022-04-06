function Zxy = mean_strength_inter(z)

[x,y] = size(z);

Zxy = (sum(sum(z)))/(x*y);
end
