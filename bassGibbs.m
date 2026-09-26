classdef bassGibbs < handle
    % BASS as one block of a larger Gibbs sampler.
    %
    % bass fits a fixed response from start to finish.  When the response is
    % itself a latent variable that another part of a sampler updates between
    % visits (data augmentation), the BASS chain has to carry its state from
    % one visit to the next and be advanced a few steps at a time against
    % whatever the response currently is.  This class does that:
    %
    %   g = bassGibbs(X, y0, nstore, 'g1', 2, 'g2', 0.01);
    %   for it = 1:n_iter
    %       ...update the latent response y...
    %       g.setResponse(y);        % keep the basis functions, swap y
    %       g.step(3);               % 3 RJMCMC moves, each followed by Gibbs
    %                                % draws of s2, beta, lam and tau
    %       mu = g.fitted();         % regression function at the training X
    %       s2 = g.s2();             % its error variance
    %       if keep, g.save(); end   % record the state as a posterior draw
    %   end
    %   mod = g.model();             % BassModel: predict(), plot()
    %
    % Validity.  Each step is the ordinary BASS update: an RJ move on the
    % basis functions given (s2, tau, lam), then Gibbs draws of s2 (beta
    % integrated out), beta, lam and tau.  That kernel leaves
    % p(basis, s2, beta, lam, tau | y) invariant, so interleaving it with
    % updates of y is a valid Gibbs sampler for the joint model.  Unlike
    % BassState.update, step() runs the Gibbs draws even when the RJ proposal
    % is invalid, so that beta always reflects the current response.
    %
    % Options are those of bass that describe the model or the proposals
    % (w1, w2, maxInt, maxBasis, npart, g1, g2, s2_lower, h1, h2, a_tau,
    % b_tau) plus the options added for this use:
    %   intercept     include a constant basis function (default true)
    %   center_basis  center basis functions over the training inputs
    %                 (default false).  intercept = false with
    %                 center_basis = true gives a regression function that
    %                 averages to exactly zero over the training inputs.
    %   s2_fixed      hold s2 fixed at this value (default NaN: sample it)
    % nstore is the number of states save() can record.

    properties
        bm      % BassModel: data, prior, the chain's state and saved samples
    end

    methods
        function obj = bassGibbs(xx, y, nstore, options)
            arguments
                xx {mustBeNumeric}
                y {mustBeNumeric}
                nstore (1,1) {mustBePositive, mustBeInteger} = 1
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
                options.intercept = true
                options.center_basis = false
                options.s2_fixed = NaN
            end
            y = y(:);
            b_tau = options.b_tau;
            if isnan(b_tau), b_tau = length(y) / 2; end
            npart = options.npart;
            if isnan(npart), npart = min(20, .1 * length(y)); end
            bd = BassData(xx, y);
            maxInt = min(options.maxInt, bd.p);
            bp = BassPrior(maxInt, options.maxBasis, npart, options.g1, options.g2, ...
                options.s2_lower, options.h1, options.h2, options.a_tau, b_tau, ...
                options.w1, options.w2);
            bp.intercept = options.intercept;
            bp.center_basis = options.center_basis;
            bp.s2_fixed = options.s2_fixed;
            obj.bm = BassModel(bd, bp, nstore);
            obj.bm.state.gibbsStep();       % so beta exists before the first step
        end

        function setResponse(obj, y)
            % Replace the response, keeping the current basis functions.
            y = y(:);
            if numel(y) ~= obj.bm.data.n
                error('bassGibbs:size', 'y must have %d elements.', obj.bm.data.n);
            end
            obj.bm.state.setResponse(y);
            obj.bm.data.y = y;
            obj.bm.data.ssy = sum(y .^ 2);
        end

        function step(obj, n_steps)
            % Advance the chain by n_steps updates against the current response.
            if nargin < 2, n_steps = 1; end
            for k = 1:n_steps
                obj.bm.state.rjStep();
                obj.bm.state.gibbsStep();
            end
        end

        function mu = fitted(obj)
            % Current regression function at the training inputs (n x 1).
            mu = obj.bm.state.fitted();
        end

        function s2 = s2(obj)
            % Current error variance.
            s2 = obj.bm.state.s2;
        end

        function n = nbasis(obj)
            n = obj.bm.state.nbasis;
        end

        function save(obj)
            % Record the current state as the next posterior draw.
            if obj.bm.k > obj.bm.nstore
                error('bassGibbs:full', 'All %d draws have been saved.', obj.bm.nstore);
            end
            obj.bm.writeState();
        end

        function mod = model(obj)
            % The saved draws as a BassModel, for predict() and plot().
            n_saved = obj.bm.k - 1;
            if n_saved < obj.bm.nstore
                warning('bassGibbs:partial', ...
                    'Only %d of %d draws saved; predict with mcmc_use = 1:%d.', ...
                    n_saved, obj.bm.nstore, n_saved);
            end
            mod = obj.bm;
        end
    end
end
