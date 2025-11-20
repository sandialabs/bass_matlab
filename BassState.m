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
            obj.s2 = 1;
            obj.nbasis = 0;
            obj.tau = 1;
            obj.s2_rate = 1;
            obj.R = 1;
            obj.lam = 1;
            obj.I_star = ones(prior.maxInt,1) * prior.w1;
            obj.I_vec = obj.I_star / sum(obj.I_star);
            obj.z_star = ones(data.p,1) .* prior.w2;
            obj.z_vec = obj.z_star ./ sum(obj.z_star);
            obj.basis = ones(data.n, 1);
            obj.nc = 1;
            obj.knots = zeros(prior.maxBasis, prior.maxInt);
            obj.knots_ind = zeros(prior.maxBasis, prior.maxInt);
            obj.signs = zeros(prior.maxBasis, prior.maxInt);
            obj.vs = zeros(prior.maxBasis, prior.maxInt);
            obj.n_int = zeros(prior.maxBasis, 1);
            obj.Xty = zeros(prior.maxBasis + 2, 1);
            obj.Xty(1) = sum(data.y);
            obj.XtX = zeros(prior.maxBasis + 2, prior.maxBasis + 2);
            obj.XtX(1, 1) = data.n;
            obj.R = [sqrt(data.n)];
            obj.R_inv_t = 1 ./ sqrt(data.n);
            obj.bhat = mean(data.y);
            obj.qf = (sqrt(data.n) * mean(data.y)).^2;
            obj.count = zeros(3,1);
            obj.cmod = false;  % has the state changed since the last write (i.e., has a birth, death, or change been accepted)?
        end

        function obj = log_post(obj)
            % get current log posterior

            lp1 = (- (obj.s2_rate + obj.prior.g2) / self.s2 ...
                - (obj.data.n / 2 + 1 + (obj.nbasis + 1) / 2 + obj.prior.g1) * log(obj.s2) ...
                + sum(log(abs(diag(obj.R)))) ...
                + (obj.prior.a_tau + (obj.nbasis + 1) / 2 - 1) * log(obj.tau) - obj.prior.a_tau * obj.tau ...
                - (obj.nbasis + 1) / 2 * log(2 * pi) ...
                + (obj.prior.h1 + obj.nbasis - 1) * log(obj.lam) - obj.lam * (obj.prior.h2 + 1));

            obj.lp = lp1;
        end

        function obj = update(obj)
            % Update the current state using a RJMCMC step
            % (and Gibbs steps at the end of this function)

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

                alpha = .5 / obj.s2 * (qf_cand.qf - obj.qf) / (1 + obj.tau) + log(obj.lam) - log(obj.nc) ...
                    + log(1 / 3) - log(1 / 3) - cand.lbmcmp + .5 * log(obj.tau) - .5 * log(1 + obj.tau);

                if log(rand) < alpha
                    obj.cmod = true;
                    % note, XtX and Xty are already updated
                    obj.nbasis = obj.nbasis + 1;
                    obj.nc = obj.nbasis + 1;
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
                ind(tokill_ind+1) = [];

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
                    obj.nc = obj.nbasis + 1;
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

                    obj.basis(:,tokill_ind+1) = [];
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

                ata = cand.basis' * cand.basis;
                Xta = obj.basis' * cand.basis;
                aty = cand.basis' * obj.data.y;

                ind = 1:obj.nc;
                XtX_cand = obj.XtX(ind,ind);
                XtX_cand(tochange_basis+1, :) = Xta;
                XtX_cand(:, tochange_basis+1) = Xta;
                XtX_cand(tochange_basis+1, tochange_basis+1) = ata;

                Xty_cand = obj.Xty(1:obj.nc);
                Xty_cand(tochange_basis+1) = aty;

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

                    obj.basis(:, tochange_basis+1) = cand.basis;
                end
            end

            a_s2 = obj.prior.g1 + obj.data.n / 2;
            b_s2 = obj.prior.g2 + .5 * (obj.data.ssy - (obj.bhat' * obj.Xty(1:obj.nc)) / (1 + obj.tau));
            if b_s2 < 0
                obj.prior.g2 = obj.prior.g2 + 1.e-10;
                b_s2 = obj.prior.g2 + .5 * (obj.data.ssy - (obj.bhat' * obj.Xty(1:obj.nc)) / (1 + obj.tau));
            end
            obj.s2 = 1 / gamrnd(a_s2, 1 / b_s2);

            obj.beta = obj.bhat / (1 + obj.tau) + (obj.R_inv_t * randn(obj.nc,1)) * sqrt(obj.s2 / (1 + obj.tau));

            a_lam = obj.prior.h1 + obj.nbasis;
            b_lam = obj.prior.h2 + 1;
            obj.lam = gamrnd(a_lam, 1 / b_lam);

            temp = obj.R * obj.beta;
            qf2 = temp' * temp;
            a_tau = obj.prior.a_tau + (obj.nbasis + 1) / 2;
            b_tau = obj.prior.b_tau + .5 * qf2 / obj.s2;
            obj.tau = gamrnd(a_tau, 1 / b_tau);
        end
    end
end

