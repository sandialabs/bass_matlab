classdef BassBasis
    % Structure for functional response BASS model using a basis decomposition, gets a list of BASS models

    properties
        basis
        xx
        y
        newy
        y_mean
        y_sd
        trunc_error
        nbasis
        bm_list
    end

    methods
        function obj = BassBasis(xx, y, basis, newy, y_mean, y_sd, trunc_error, ncores, opts, fit)
            % Fit BASS model with multivariate/functional response by projecting onto user specified basis.

            % xx: matrix  of predictors of dimension nxp, where n is the number of training examples and
            %     p is the number of inputs (features).
            % y: response matrix of dimension nxq, where q is the number of multivariate/functional
            %    responses.
            % basis: matrix of basis functions of dimension qxk.
            % newy: matrix of y projected onto basis, dimension kxn.
            % y_mean: vector of length q with the mean if y was centered before obtaining newy.
            % y_sd: vector of length q with the standard deviation if y was scaled before obtaining newy.
            % trunc_error: numpy array of projection truncation errors (dimension qxn)
            % ncores: number of threads to use when fitting independent BASS models (integer less than or equal to
            %         npc).
            % opts: optional arguments to bass function.

            arguments
                xx
                y
                basis
                newy = nan
                y_mean = nan
                y_sd = nan
                trunc_error = nan
                ncores = nan
                opts = nan
                fit = true
            end

            obj.basis = basis;
            obj.xx = xx;
            obj.y = y;
            obj.newy = newy;
            obj.y_mean = y_mean;
            obj.y_sd = y_sd;
            obj.trunc_error = trunc_error;
            obj.nbasis = size(basis,2);

            if fit
                if ncores == 1
                    obj.bm_list = cell(1,obj.nbasis);
                    for ii = 1:obj.nbasis
                        obj.bm_list{ii} = bass(obj.xx, obj.newy(ii,:)', opts.nmcmc, opts.nburn, opts.thin, ...
                            opts.w1, opts.w2, opts.maxInt, opts.maxBasis, opts.npart, opts.g1, ...
                            opts.g2, opts.s2_lower, opts.h1, opts.h2, opts.a_tau, opts.b_tau, opts.verbose);
                    end
                else
                    temp = PoolBass(obj.xx, obj.newy, opts);
                    obj.bm_list = temp.fit(ncores, obj.nbasis);
                end
            else
                obj.bm_list = cell(1,obj.nbasis);
            end
        end

        function out2 = predict(obj, X, mcmc_use, nugget, trunc_error, ncores)
            % Predict the functional response at new inputs.

            % X: matrix of predictors with dimension nxp, where n is the number of prediction points and
            %    p is the number of inputs (features). p must match the number of training inputs, and the order of the
            %    columns must also match.
            % mcmc_use: which MCMC samples to use (list of integers of length m).  Defaults to all MCMC samples.
            % nugget: whether to use the error variance when predicting.  If False, predictions are for mean function.
            % trunc_error: whether to use truncation error when predicting.
            % ncores: number of cores to use while predicting (integer).  In almost all cases, use ncores=1.
            % returns a numpy array of predictions with dimension mxnxq, with first dimension corresponding to MCMC samples,
            % second dimension corresponding to prediction points, and third dimension corresponding to
            % multivariate/functional response.

            arguments
                obj BassBasis
                X {mustBeNumeric}
                mcmc_use = NaN
                nugget = false
                trunc_error = false
                ncores = 1
            end

            if isnan(mcmc_use)
                mcmc_use = 1:obj.bm_list{1}.nstore;
            end

            if ncores == 1
                pred_coefs = cell(1,obj.nbasis);
                for i = 1:obj.nbasis
                    pred_coefs{i} = obj.bm_list{i}.predict(X, mcmc_use, nugget);
                end
            else
                temp = PoolBassPredict(X, mcmc_use, nugget, obj.bm_list);
                pred_coefs = temp.predict(ncores, obj.nbasis);
            end

            out = zeros(size(pred_coefs{1},2),size(obj.basis,1),length(mcmc_use));
            tmp = cat(3,pred_coefs{:});
            for i = 1:length(mcmc_use)
                tmp1 = squeeze(tmp(i,:,:));
                if size(tmp,2) == 1
                    tmp1 = tmp1';
                end
                out(:,:,i) = tmp1 * obj.basis';
            end
            tmp_mean = repmat(obj.y_mean,size(pred_coefs{1},2),1,length(mcmc_use));
            tmp_sd = repmat(obj.y_sd,size(pred_coefs{1},2),1,length(mcmc_use));
            out2 = out .* tmp_sd + tmp_mean;
            if trunc_error
                out2 = out2 + reshape(obj.trunc_error(:, randsample(1:size(obj.trunc_error,2), size(out,1)*size(out,3), true)), size(out));
            end

        end

        function plot(obj)
            % Trace plots and predictions/residuals

            % * top left - trace plot of number of basis functions (excluding burn-in and thinning)
            % * top right - trace plot of residual variance
            % * bottom left - training data against predictions
            % * bottom right - histogram of residuals (posterior mean) with assumed Gaussian overlaid.

            figure()
            subplot(2,2,1)
            hold all
            for i = 1:obj.nbasis
                plot(obj.bm_list{i}.samples.nbasis)
            end
            ylabel('number of basis functions')
            xlabel('MCMC iteration (post-burn)')

            subplot(2,2,2)
            hold all
            for i = 1:obj.nbasis
                plot(obj.bm_list{i}.samples.s2)
            end
            ylabel('error variance')
            xlabel('MCMC iteration (post-burn)')

            subplot(2,2,3)
            yhat = mean(obj.predict(obj.xx),3);
            scatter(obj.y, yhat)
            refline(1,0)
            xlabel('observed')
            ylabel('posterior prediction')

            subplot(2,2,4)
            tmp = obj.y-yhat;
            histfit(reshape(tmp,numel(yhat),1))
            xlabel('residuals')
            ylabel('density')
        end
    end
end
