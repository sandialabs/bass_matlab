classdef BassState < handle
    % The current state of the RJMCMC chain, with methods for getting the
    % log posterior and for updating the state

    properties
        data
        prior
        s2
        nbasis
        tau
        s2_rate
        R
        lam
        I_star
        I_vec
        z_star
        z_vec
        basis
        nc
        knots
        knots_ind
        signs
        vs
        n_int
        Xty
        XtX
        R_inv_t
        bhat
        qf
        count
        cmod
        lp
        beta
    end

    methods
        function obj = BassState(data, prior)
            obj.data = data;
            obj.prior = prior;
            if isnan(prior.s2_fixed)
                obj.s2 = 1;
            else
                obj.s2 = prior.s2_fixed;
            end
            obj.nbasis = 0;
            obj.tau = 1;
            obj.s2_rate = 1;
            obj.R = 1;
            obj.lam = 1;
            obj.I_star = ones(prior.maxInt,1) * prior.w1;
            obj.I_vec = obj.I_star / sum(obj.I_star);
            obj.z_star = ones(data.p,1) .* prior.w2;
            obj.z_vec = obj.z_star ./ sum(obj.z_star);
            % Column 1 of the basis is the intercept, absent when
            % prior.intercept is false.
            i0 = double(prior.intercept);
            obj.basis = ones(data.n, i0);
            obj.nc = i0;
            obj.knots = zeros(prior.maxBasis, prior.maxInt);
            obj.knots_ind = zeros(prior.maxBasis, prior.maxInt);
            obj.signs = zeros(prior.maxBasis, prior.maxInt);
            obj.vs = zeros(prior.maxBasis, prior.maxInt);
            obj.n_int = zeros(prior.maxBasis, 1);
            obj.Xty = zeros(prior.maxBasis + 2, 1);
            obj.XtX = zeros(prior.maxBasis + 2, prior.maxBasis + 2);
            if prior.intercept
                obj.Xty(1) = sum(data.y);
                obj.XtX(1, 1) = data.n;
                obj.R = [sqrt(data.n)];
                obj.R_inv_t = 1 ./ sqrt(data.n);
                obj.bhat = mean(data.y);
                obj.qf = (sqrt(data.n) * mean(data.y)).^2;
            else
                obj.R = zeros(0, 0);
                obj.R_inv_t = zeros(0, 0);
                obj.bhat = zeros(0, 1);
                obj.qf = 0;
            end
            obj.count = zeros(3,1);
            % Has the basis changed since the last write?  True at the start so
            % that the initial (empty) basis is recorded if no move has been
            % accepted before the first write.
            obj.cmod = true;
        end

        function obj = log_post(obj)
            % get current log posterior

            lp1 = (- (obj.s2_rate + obj.prior.g2) / obj.s2 ...
                - (obj.data.n / 2 + 1 + obj.nc / 2 + obj.prior.g1) * log(obj.s2) ...
                + sum(log(abs(diag(obj.R)))) ...
                + (obj.prior.a_tau + obj.nc / 2 - 1) * log(obj.tau) - obj.prior.a_tau * obj.tau ...
                - obj.nc / 2 * log(2 * pi) ...
                + (obj.prior.h1 + obj.nbasis - 1) * log(obj.lam) - obj.lam * (obj.prior.h2 + 1));

            obj.lp = lp1;
        end

        function obj = update(obj)
            % Update the current state using a RJMCMC step, then Gibbs steps
            % for s2, beta, lam and tau.  As before, the Gibbs steps are
            % skipped when the RJ proposal is invalid (too few nonzero
            % points, or a rank-deficient basis).
            if obj.rjStep()
                obj.gibbsStep();
            end
        end

        function ok = rjStep(obj)
            % One RJMCMC birth/death/change move.  Returns false if the
            % proposal was invalid and nothing was evaluated.
            ok = false;
            i0 = double(obj.prior.intercept);
            move_type = randsample(1:3,1);

            if obj.nbasis == 0
                move_type = 1;
            end

            if obj.nbasis == obj.prior.maxBasis
                move_type = randsample(2:3, 1);
            end

            if move_type == 1
                % BIRTH step

                cand = genCandBasis(obj.prior.maxInt, obj.I_vec, obj.z_vec, obj.data.p, obj.data.xx);

                if sum(cand.basis(:) > 0) < obj.prior.npart
                    return
                end
                if obj.prior.center_basis
                    cand.basis = cand.basis - mean(cand.basis);
                end

                ata = cand.basis' * cand.basis;
                Xta = obj.basis' * cand.basis;
                aty = cand.basis' * obj.data.y;

                obj.Xty(obj.nc+1) = aty;
                obj.XtX(1:(obj.nc), obj.nc+1) = Xta;
                obj.XtX(obj.nc+1, 1:(obj.nc)) = Xta;
                obj.XtX(obj.nc+1, obj.nc+1) = ata;

                qf_cand = getQf(obj.XtX(1:(obj.nc+1), 1:(obj.nc+1)), obj.Xty(1:(obj.nc+1)));

                if ~qf_cand.fullrank
                    return
                end

                alpha = .5 / obj.s2 * (qf_cand.qf - obj.qf) / (1 + obj.tau) + log(obj.lam) - log(obj.nbasis + 1) ...
                    + log(1 / 3) - log(1 / 3) - cand.lbmcmp + .5 * log(obj.tau) - .5 * log(1 + obj.tau);

                if log(rand) < alpha
                    obj.cmod = true;
                    % note, XtX and Xty are already updated
                    obj.nbasis = obj.nbasis + 1;
                    obj.nc = obj.nbasis + i0;
                    obj.qf = qf_cand.qf;
                    obj.bhat = qf_cand.bhat;
                    obj.R = qf_cand.R;
                    obj.R_inv_t = obj.R\eye(obj.nc);
                    obj.count(1) = obj.count(1) + 1;
                    obj.n_int(obj.nbasis) = cand.n_int;
                    obj.knots(obj.nbasis, 1:(cand.n_int)) = cand.knots;
                    obj.knots_ind(obj.nbasis, 1:(cand.n_int)) = cand.knots_ind;
                    obj.signs(obj.nbasis, 1:(cand.n_int)) = cand.signs;
                    obj.vs(obj.nbasis, 1:(cand.n_int)) = cand.vs;

                    obj.I_star(cand.n_int) = obj.I_star(cand.n_int) + 1;
                    obj.I_vec = obj.I_star / sum(obj.I_star);
                    obj.z_star(cand.vs) = obj.z_star(cand.vs) + 1;
                    obj.z_vec = obj.z_star / sum(obj.z_star);

                    obj.basis = [obj.basis, cand.basis];
                end

            elseif move_type == 2
                % DEATH step

                tokill_ind = randsample(obj.nbasis,1);
                ind = 1:obj.nc;
                ind(tokill_ind+i0) = [];

                qf_cand = getQf(obj.XtX(ind,ind), obj.Xty(ind));

                if ~qf_cand.fullrank
                    return
                end

                I_star1 = obj.I_star;
                I_star1(obj.n_int(tokill_ind)) = I_star1(obj.n_int(tokill_ind)) - 1;
                I_vec1 = I_star1 / sum(I_star1);
                z_star1 = obj.z_star;
                z_star1(obj.vs(tokill_ind, 1:obj.n_int(tokill_ind))) = z_star1(obj.vs(tokill_ind,1:obj.n_int(tokill_ind))) - 1;

                z_vec1 = z_star1 / sum(z_star1);

                lbmcmp = logProbChangeMod(obj.n_int(tokill_ind), obj.vs(tokill_ind, 1:obj.n_int(tokill_ind)), I_vec1, ...
                    z_vec1, obj.data.p, obj.prior.maxInt);

                alpha = .5 / obj.s2 * (qf_cand.qf - obj.qf) / (1 + obj.tau) - log(obj.lam) + log(obj.nbasis) ...
                    + log(1 / 3) - log(1 / 3) + lbmcmp - .5 * log(obj.tau) + .5 * log(1 + obj.tau);

                if log(rand) < alpha
                    obj.cmod = true;
                    obj.nbasis = obj.nbasis - 1;
                    obj.nc = obj.nbasis + i0;
                    obj.qf = qf_cand.qf;
                    obj.bhat = qf_cand.bhat;
                    obj.R = qf_cand.R;
                    obj.R_inv_t = obj.R\eye(obj.nc);
                    obj.count(2) = obj.count(2) + 1;

                    obj.Xty(1:obj.nc) = obj.Xty(ind);
                    obj.XtX(1:obj.nc, 1:obj.nc) = obj.XtX(ind, ind);

                    temp = obj.n_int(1:(obj.nbasis+1));
                    temp(tokill_ind) = [];
                    obj.n_int = obj.n_int * 0;
                    obj.n_int(1:(obj.nbasis)) = temp(:);

                    temp = obj.knots(1:(obj.nbasis+1), :);
                    temp(tokill_ind,:) = [];
                    obj.knots = obj.knots * 0;
                    obj.knots(1:(obj.nbasis), :) = temp;

                    temp = obj.knots_ind(1:(obj.nbasis+1), :);
                    temp(tokill_ind,:) = [];
                    obj.knots_ind = obj.knots_ind * 0;
                    obj.knots_ind(1:(obj.nbasis), :) = temp;

                    temp = obj.signs(1:(obj.nbasis+1), :);
                    temp(tokill_ind,:) = [];
                    obj.signs = obj.signs * 0;
                    obj.signs(1:(obj.nbasis), :) = temp;

                    temp = obj.vs(1:(obj.nbasis + 1), :);
                    temp(tokill_ind,:) = [];
                    obj.vs = obj.vs * 0;
                    obj.vs(1:(obj.nbasis), :) = temp;

                    obj.I_star = I_star1;
                    obj.I_vec = I_vec1;
                    obj.z_star = z_star1;
                    obj.z_vec = z_vec1;

                    obj.basis(:,tokill_ind+i0) = [];
                end

            else
                % CHANGE step
                tochange_basis = randsample(obj.nbasis,1);
                tochange_int = randsample(obj.n_int(tochange_basis),1);

                cand = genBasisChange(obj.knots(tochange_basis, 1:obj.n_int(tochange_basis)), ...
                    obj.signs(tochange_basis, 1:obj.n_int(tochange_basis)), ...
                    obj.vs(tochange_basis, 1:obj.n_int(tochange_basis)), ...
                    obj.knots_ind(tochange_basis, 1:obj.n_int(tochange_basis)), tochange_int, obj.data.xx);

                if sum(cand.basis > 0) < obj.prior.npart
                    return
                end
                if obj.prior.center_basis
                    cand.basis = cand.basis - mean(cand.basis);
                end

                ata = cand.basis' * cand.basis;
                Xta = obj.basis' * cand.basis;
                aty = cand.basis' * obj.data.y;

                ind = 1:obj.nc;
                XtX_cand = obj.XtX(ind,ind);
                XtX_cand(tochange_basis+i0, :) = Xta;
                XtX_cand(:, tochange_basis+i0) = Xta;
                XtX_cand(tochange_basis+i0, tochange_basis+i0) = ata;

                Xty_cand = obj.Xty(1:obj.nc);
                Xty_cand(tochange_basis+i0) = aty;

                qf_cand = getQf(XtX_cand, Xty_cand);

                if ~qf_cand.fullrank
                    return
                end

                alpha = .5 / obj.s2 * (qf_cand.qf - obj.qf) / (1 + obj.tau);

                if log(rand) < alpha
                    obj.cmod = true;
                    obj.qf = qf_cand.qf;
                    obj.bhat = qf_cand.bhat;
                    obj.R = qf_cand.R;
                    obj.R_inv_t = obj.R\eye(obj.nc);
                    obj.count(3) = obj.count(3) + 1;

                    obj.Xty(1:obj.nc) = Xty_cand;
                    obj.XtX(1:obj.nc, 1:obj.nc) = XtX_cand;

                    obj.knots(tochange_basis, 1:obj.n_int(tochange_basis)) = cand.knots;
                    obj.knots_ind(tochange_basis, 1:obj.n_int(tochange_basis)) = cand.knots_ind;
                    obj.signs(tochange_basis, 1:obj.n_int(tochange_basis)) = cand.signs;

                    obj.basis(:, tochange_basis+i0) = cand.basis;
                end
            end
            ok = true;
        end

        function obj = gibbsStep(obj)
            % Gibbs updates of s2 (beta integrated out), beta, lam and tau
            % given the current basis.
            if isnan(obj.prior.s2_fixed)
                obj = obj.sampleS2();
            else
                obj.s2 = obj.prior.s2_fixed;
            end

            obj.beta = obj.bhat / (1 + obj.tau) + (obj.R_inv_t * randn(obj.nc,1)) * sqrt(obj.s2 / (1 + obj.tau));

            a_lam = obj.prior.h1 + obj.nbasis;
            b_lam = obj.prior.h2 + 1;
            obj.lam = gamrnd(a_lam, 1 / b_lam);

            temp = obj.R * obj.beta;
            qf2 = temp' * temp;
            a_tau = obj.prior.a_tau + obj.nc / 2;
            b_tau = obj.prior.b_tau + .5 * qf2 / obj.s2;
            obj.tau = gamrnd(a_tau, 1 / b_tau);
        end

        function obj = sampleS2(obj)
            a_s2 = obj.prior.g1 + obj.data.n / 2;
            b_s2 = obj.prior.g2 + .5 * (obj.data.ssy - (obj.bhat' * obj.Xty(1:obj.nc)) / (1 + obj.tau));
            if b_s2 < 0
                obj.prior.g2 = obj.prior.g2 + 1.e-10;
                b_s2 = obj.prior.g2 + .5 * (obj.data.ssy - (obj.bhat' * obj.Xty(1:obj.nc)) / (1 + obj.tau));
            end
            obj.s2 = 1 / gamrnd(a_s2, 1 / b_s2);
        end

        function obj = setResponse(obj, y)
            % Point the chain at a new response vector, keeping the basis
            % functions.  For use inside a larger Gibbs sampler (bassGibbs),
            % where the response is a latent variable that changes between
            % updates.  Only the quantities that involve y change; beta is
            % refreshed by the next gibbsStep.
            obj.data.y = y(:);
            obj.data.ssy = sum(obj.data.y .^ 2);
            nc = obj.nc;
            obj.Xty(1:nc) = obj.basis' * obj.data.y;
            Qf = getQf(obj.XtX(1:nc, 1:nc), obj.Xty(1:nc));
            obj.bhat = Qf.bhat;
            obj.qf = Qf.qf;
        end

        function mu = fitted(obj)
            % Current regression function at the training inputs.
            mu = obj.basis * obj.beta;
        end
    end
end

