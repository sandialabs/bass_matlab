function out = logProbChangeMod(n_int, vars_use, I_vec, z_vec, p, maxInt)
% Get reversibility factor for RJMCMC acceptance ratio, and also prior

if n_int == 1
    out = (log(I_vec(n_int)) - log(2*p) + log(2*p) + log(maxInt));
else
    x = zeros(p,1);
    x(vars_use) = 1;
    lprob_vars_noReplace = log(dmwnchBass(z_vec, vars_use));
    out = (log(I_vec(n_int)) + lprob_vars_noReplace - n_int * log(2) + n_int * log(2) + log(nchoosek(p, n_int)) + log(maxInt));
end
