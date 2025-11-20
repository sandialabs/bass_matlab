function CandidateBasis = genCandBasis(maxInt, I_vec, z_vec, p, xdata)
% Generate a candidate basis for birth step, as well as the RJMCMC reversibility factor and prior

n_int = randsample(1:maxInt, 1, true, I_vec);
signs = randsample([-1, 1], n_int, true);

knots = zeros(n_int,1);
knots_ind = zeros(n_int,1);
if n_int == 1
    vs = datasample(1:p,1);
    [knots,knots_ind] = datasample(xdata(:, vs), 1);
else
    vs = sort(datasample(1:p, n_int, 'Replace', false, 'Weights', z_vec));
    for i = 1:n_int
        [knots(i), knots_ind(i)] = datasample(xdata(:,vs(i)), 1);
    end
end
basis = makeBasis(signs, vs, knots, xdata);
lbmcmp = logProbChangeMod(n_int, vs, I_vec, z_vec, p, maxInt);

CandidateBasis.basis = basis;
CandidateBasis.n_int = n_int;
CandidateBasis.signs = signs;
CandidateBasis.vs = vs;
CandidateBasis.knots = knots;
CandidateBasis.knots_ind = knots_ind;
CandidateBasis.lbmcmp = lbmcmp;
