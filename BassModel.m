classdef BassModel < handle
    % The model structure, including the current RJMCMC state and previous saved states; with methods for saving the
    % state, plotting MCMC traces, and predicting

    properties
        data
        prior
        state
        nstore
        samples
        k
        k_mod
        model_lookup
    end

    methods
        function obj = BassModel(data, prior, nstore)
            obj.data = data;
            obj.prior = prior;
            obj.state = BassState(obj.data, obj.prior);
            obj.nstore = nstore;
            s2 = zeros(nstore,1);
            lam = zeros(nstore,1);
            tau = zeros(nstore,1);
            nbasis = zeros(nstore,1);
            nbasis_models = zeros(nstore, 1);
            n_int = zeros(obj.nstore, obj.prior.maxBasis);
            signs = zeros(obj.nstore, obj.prior.maxBasis, obj.prior.maxInt);
            vs = zeros(obj.nstore, obj.prior.maxBasis, obj.prior.maxInt);
            knots = zeros(obj.nstore, obj.prior.maxBasis, obj.prior.maxInt);
            knots_ind = zeros(obj.nstore, obj.prior.maxBasis, obj.prior.maxInt);
            beta = zeros(obj.nstore, obj.prior.maxBasis + 1);
            Samples.s2 = s2;
            Samples.lam = lam;
            Samples.tau = tau;
            Samples.nbasis = nbasis;
            Samples.nbasis_models = nbasis_models;
            Samples.n_int = n_int;
            Samples.signs = signs;
            Samples.vs = vs;
            Samples.knots = knots;
            Samples.knots_ind = knots_ind;
            Samples.beta = beta;
            obj.samples = Samples;
            obj.k = 1;
            obj.k_mod = 0;
            obj.model_lookup = zeros(nstore,1);
        end

        function obj = writeState(obj)
            % Take relevant parts of state and write to storage (only manipulates storage vectors created in init)
            obj.samples.s2(obj.k) = obj.state.s2;
            obj.samples.lam(obj.k) = obj.state.lam;
            obj.samples.tau(obj.k) = obj.state.tau;
            obj.samples.beta(obj.k, 1:(obj.state.nbasis + 1)) = obj.state.beta;
            obj.samples.nbasis(obj.k) = obj.state.nbasis;

            if obj.state.cmod % basis part of state was changed
                obj.k_mod = obj.k_mod + 1;
                obj.samples.nbasis_models(obj.k_mod) = obj.state.nbasis;
                obj.samples.n_int(obj.k_mod, 1:obj.state.nbasis) = obj.state.n_int(1:obj.state.nbasis);
                obj.samples.signs(obj.k_mod, 1:obj.state.nbasis, :) = obj.state.signs(1:obj.state.nbasis, :);
                obj.samples.vs(obj.k_mod, 1:obj.state.nbasis, :) = obj.state.vs(1:obj.state.nbasis, :);
                obj.samples.knots(obj.k_mod, 1:obj.state.nbasis, :) = obj.state.knots(1:obj.state.nbasis, :);
                obj.samples.knots_ind(obj.k_mod, 1:obj.state.nbasis, :) = obj.state.knots_ind(1:obj.state.nbasis, :);
                obj.state.cmod = false;
            end

            obj.model_lookup(obj.k) = obj.k_mod;
            obj.k = obj.k + 1;
        end

        function plot(obj)
            % Trace plots and predictions/residuals

            % * top left - trace plot of number of basis functions (excluding burn-in and thinning)
            % * top right - trace plot of residual variance
            % * bottom left - training data against predictions
            % * bottom right - histogram of residuals (posterior mean) with assumed Gaussian overlaid.

            figure()
            subplot(2,2,1)
            plot(obj.samples.nbasis)
            ylabel('number of basis functions')
            xlabel('MCMC iteration (post-burn)')

            subplot(2,2,2)
            plot(obj.samples.s2)
            ylabel('error variance')
            xlabel('MCMC iteration (post-burn)')

            subplot(2,2,3)
            yhat = mean(obj.predict(obj.data.xx_orig, NaN, false),1);
            scatter(obj.data.y, yhat)
            refline(1,0)
            xlabel('observed')
            ylabel('posterior prediction')

            subplot(2,2,4)
            histfit(obj.data.y(:)-yhat(:))
            xlabel('residuals')
            ylabel('density')
        end

        function mat = makeBasisMatrix(obj, model_ind, X)
            % Make basis matrix for model
            nb = obj.samples.nbasis_models(model_ind);
            n = size(X,1);
            mat = zeros(n,nb+1);
            mat(:,1) = ones(n,1);
            for m = 1:nb
                ind = 1:obj.samples.n_int(model_ind,m);
                signs = squeeze(obj.samples.signs(model_ind,m,ind));
                vs = squeeze(obj.samples.vs(model_ind,m,ind));
                knots = squeeze(obj.samples.knots(model_ind,m,ind));
                mat(:, m+1) = makeBasis(signs, vs, knots, X);
            end
        end

        function out = predict(obj, X, mcmc_use, nugget)
            % BASS prediction using new inputs (after training).

            % X: matrix of predictors with dimension nxp, where n is the number of prediction points and
            %    p is the number of inputs (features). p must match the number of training inputs, and the order of the
            %    columns must also match.
            % mcmc_use: which MCMC samples to use (vector of integers of length m).  Defaults to all MCMC samples.
            % nugget: whether to use the error variance when predicting.  If False, predictions are for mean function.
            %         a matrix of predictions with dimension mxn, with rows corresponding to MCMC samples and
            %         columns corresponding to prediction points.
            if nargin < 3
                mcmc_use = NaN;
                nugget = false;
            end

            if nargin < 4
                nugget = false;
            end

            Xs = normalizebass(X, obj.data.bounds);
            if isnan(mcmc_use)
                mcmc_use = 1:obj.nstore;
            end

            out = zeros(length(mcmc_use), size(Xs,1));
            models = obj.model_lookup(mcmc_use);
            umodels = unique(models);
            k1 = 1;
            for j = 1:length(umodels)
                mcmc_use_j = mcmc_use(models==umodels(j));
                nn = length(mcmc_use_j);
                out(k1:(k1+nn-1),:) = obj.samples.beta(mcmc_use_j, 1:(obj.samples.nbasis_models(umodels(j))+1)) * obj.makeBasisMatrix(umodels(j),Xs)';
                k1 = k1 + nn;
            end

            if nugget
                out = out + bsxfun(@times,sqrt(obj.samples.s2(mcmc_use)),randn(size(out)));
            end

        end

    end
end
