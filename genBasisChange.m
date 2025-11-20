function BasisChange = genBasisChange(knots, signs, vs, knots_ind, tochange_int, xdata)
% Generate a condidate basis for change step
knots_cand = knots;
signs_cand = signs;
knots_ind_cand = knots_ind;
signs_cand(tochange_int) = randsample([-1, 1], 1);
[knots_cand(tochange_int), knots_ind_cand(tochange_int)] = datasample(xdata(:, vs(tochange_int)), 1);

basis = makeBasis(signs_cand, vs, knots_cand, xdata);
BasisChange.basis = basis;
BasisChange.signs = signs_cand;
BasisChange.vs = vs;
BasisChange.knots = knots_cand;
BasisChange.knots_ind = knots_ind_cand;
