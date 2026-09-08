function ss = dmwnchBass(z_vec, vars_use)
% Multivariate Walenius' noncentral hypergeometric density function with some variables fixed"

z_rm = z_vec;
z_rm(vars_use) = [];
alpha = z_vec(vars_use) ./ sum(z_rm);
j = length(alpha);
ss = 1 + (-1)^j * 1 / (sum(alpha) + 1);
% Inclusion-exclusion over subsets of size 1 ... j-1; the size-j term is the
% leading (-1)^j term above.  Subset of size (i+1) carries sign (-1)^(i+1).
for i = 0:(j-2)
    idx = nchoosek(1:j,i+1);
    temp = alpha(idx);
    ss = ss + (-1)^(i+1) * sum(1./(sum(temp,2)+ 1));
end
