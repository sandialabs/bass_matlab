function bm = bassPCA(xx, y, npc, percVar, ncores, center, scale, opts)
% Wrapper to get principal components and call BassBasis, which then calls bass function to fit the BASS model for
% functional (or multivariate) response data.

% xx: matrix of predictors of dimension nxp, where n is the number of training examples and p is
%     the number of inputs (features).
% y: response matrix of dimension nxq, where q is the number of multivariate/functional
%    responses.
% npc: number of principal components to use (integer, optional if percVar is specified).
% percVar: percent (between 0 and 100) of variation to explain when choosing number of principal components
%          (if npc=None).
% ncores: number of threads to use when fitting independent BASS models (integer less than or equal to npc).
% center: whether to center the responses before principal component decomposition (boolean).
% scale: whether to scale the responses before principal component decomposition (boolean).
% opts: optional arguments to bass function.
% returns object of class BassBasis, with predict and plot functions.

arguments
    xx {mustBeNumeric}
    y {mustBeNumeric}
    npc = NaN;
    percVar = 99.9;
    ncores = 1;
    center = true;
    scale = false;
    opts.nmcmc = 10000;
    opts.nburn = 9000;
    opts.thin = 1;
    opts.w1 = 5;
    opts.w2 = 5;
    opts.maxInt = 3;
    opts.maxBasis = 1000;
    opts.npart = NaN;
    opts.g1 = 0;
    opts.g2 = 0;
    opts.s2_lower = 0;
    opts.h1 = 10;
    opts.h2 = 10;
    opts.a_tau = 0.5;
    opts.b_tau = NaN;
    opts.verbose = true;
end

setup = BassPCAsetup(y, center, scale);

if isnan(npc)
    cs = cumsum(setup.evals) / sum(setup.evals) * 100;
    npc = find(cs >= percVar, 1);
end

if ncores > npc
    ncores = npc;
end

basis = setup.basis(:, 1:npc);
newy = setup.newy(1:npc, :);
trunc_error = basis * newy - setup.y_scale';

fprintf('Starting bassPCA with %d components, using %d cores.\n',npc, ncores)

bm = BassBasis(xx, y, basis, newy, setup.y_mean, setup.y_sd, trunc_error, ncores, opts);
