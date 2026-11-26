function BasisChange = genBasisChangeLocal(knots, signs, vs, knots_ind, tochange_int, data, width)
% Candidate basis for a local change move: the knot of factor tochange_int
% moves to the data point k places away in the sorted order of its
% variable, k uniform on +-(1:width); the sign is kept.  The proposal is
% symmetric, and the knot prior is uniform over the data points, so the
% acceptance ratio is the likelihood ratio alone.  ok is false when the
% move would leave the data (a rejection, which keeps it symmetric).
v = vs(tochange_int);
r = data.rnk(knots_ind(tochange_int), v) + randi(width) * (2 * randi(2) - 3);
BasisChange.ok = r >= 1 && r <= data.n;
if ~BasisChange.ok
    return
end
knots_ind_cand = knots_ind;
knots_cand = knots;
knots_ind_cand(tochange_int) = data.ord(r, v);
knots_cand(tochange_int) = data.xx(knots_ind_cand(tochange_int), v);
BasisChange.basis = makeBasis(signs, vs, knots_cand, data.xx);
BasisChange.signs = signs;
BasisChange.vs = vs;
BasisChange.knots = knots_cand;
BasisChange.knots_ind = knots_ind_cand;
end
