function bm = bass(xx, y, options)
% **Bayesian Adaptive Spline Surfaces - model fitting**

% This function takes training data, priors, and algorithmic constants and fits a BASS model.  The result is a set of
%posterior samples of the model.  The resulting object has a predict function to generate posterior predictive
%samples.  Default settings of priors and algorithmic parameters should only be changed by users who understand
%the model.

% xx: matrix of predictors of dimension nxp, where n is the number of training examples and p is
%     the number of inputs (features).
% y: response vector of length n.
% options structure with fields
% nmcmc: total number of MCMC iterations (integer)
% nburn: number of MCMC iterations to throw away as burn-in (integer, less than nmcmc).
% thin: number of MCMC iterations to thin (integer).
% w1: nominal weight for degree of interaction, used in generating candidate basis functions. Should be greater
%     than 0.
% w2: nominal weight for variables, used in generating candidate basis functions. Should be greater than 0.
% maxInt: maximum degree of interaction for spline basis functions (integer, less than p)
% maxBasis: maximum number of tensor product spline basis functions (integer)
% npart: minimum number of non-zero points in a basis function. If the response is functional, this refers only
%        to the portion of the basis function coming from the non-functional predictors. Defaults to 20 or 0.1 times the
%        number of observations, whichever is smaller.
% g1: shape for IG prior on residual variance.
% g2: scale for IG prior on residual variance.
% s2_lower: lower bound for residual variance.
% h1: shape for gamma prior on mean number of basis functions.
% h2: scale for gamma prior on mean number of basis functions.
% a_tau: shape for gamma prior on 1/g in g-prior.
% b_tau: scale for gamma prior on 1/g in g-prior.
% verbose: boolean for printing progress
% returns an object of class BassModel, which includes predict and plot functions.

arguments
    xx {mustBeNumeric}
    y {mustBeNumeric}
    options.nmcmc = 10000
    options.nburn = 9000
    options.thin = 1
    options.w1 = 5
    options.w2 = 5
    options.maxInt = 3
    options.maxBasis = 1000
    options.npart = NaN
    options.g1 = 0
    options.g2 = 0
    options.s2_lower = 0
    options.h1 = 10
    options.h2 = 10
    options.a_tau = 0.5
    options.b_tau = NaN
    options.verbose = true
end

nmcmc = options.nmcmc;
nburn = options.nburn;
thin = options.thin;
w1 = options.w1;
w2 = options.w2;
maxInt = options.maxInt;
maxBasis = options.maxBasis;
npart = options.npart;
g1 = options.g1;
g2 = options.g2;
s2_lower = options.s2_lower;
h1 = options.h1;
h2 = options.h2;
a_tau = options.a_tau;
b_tau = options.b_tau;
verbose = options.verbose;

if isnan(b_tau)
    b_tau = length(y)/2;
end
if isnan(npart)
    npart = min(20,.1*length(y));
end
bd = BassData(xx,y);
if bd.p < maxInt
    maxInt = bd.p;
end
if verbose
    obj = ProgressBar(nmcmc, 'Title', 'Running BASS MCMC');
end
bp = BassPrior(maxInt, maxBasis, npart, g1, g2, s2_lower, h1, h2, a_tau, b_tau, w1, w2);
nstore = floor((nmcmc-nburn) / thin);
bm = BassModel(bd, bp, nstore);
for i = 1:nmcmc
    bm.state.update();
    if i > (nburn - 1) && mod(i - nburn, thin) == 0
        bm.writeState();
    end
    if verbose
        obj.step([], [], []);
    end
end
if verbose
    obj.release()
end
