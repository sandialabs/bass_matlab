classdef sobolBasis
    % BASS Sensitivity Analysis

    properties
        mod
        S
        S_var
        T_var
        Var_tot
        names_ind
        xx
    end

    methods
        function obj = sobolBasis(mod)
            arguments
                mod BassBasis
            end
            obj.mod = mod;
        end

        function obj = decomp(obj, int_order, prior, mcmc_use, nind, ncores)
            % Decomposes the variance of the BASS model into variance due to main effects,
            % two way interactions, and so on, similar to the ANOVA decomposition for linear models.
            % Uses the Sobol' decomposition, which can be done analytically for MARS models.
            %
            % Parameters
            % int_order: an integer indicating the highest order of interactions to include in the Sobol decomposition.
            % prior: a cell array with the same number of elements as there are inputs to mod.
            %        Each element specifies the prior for the particular input.  Each prior is specified as a
            %        struct with elements (one of "normal", "student", or "uniform"), "trunc" (a vector of dimension 2
            %        indicating the lower and upper truncation bounds, taken to be the data bounds if omitted), and for "normal"
            %        or "student" priors, "mean" (scalar mean of the Normal/Student, or a vector of means for a mixture of
            %        Normals or Students), "sd" (scalar standard deviation of the Normal/Student, or a vector of standard
            %        deviations for a mixture of Normals or Students), "df" (scalar degrees of freedom of the Student,
            %        or a vector of degrees of freedom for a mixture of Students), and "weights" (a vector of weights that
            %        sum to one for the mixture components, or the scalar 1).  If unspecified, a uniform is assumed with the same
            %        bounds as are represented in the input to mod.
            % mcmc_use: an integer indicating which MCMC iteration to use for sensitivity analysis. Defaults to the last iteration.
            % nind: number of Sobol indices to keep (will keep the largest nind).
            % ncores: number of cores to use (default = 1)

            arguments
                obj
                int_order
                prior = NaN
                mcmc_use = NaN
                nind = NaN
                ncores = 1
            end

            if ncores > 1
                if isempty(gcp('nocreate'))
                    parpool(ncores);
                end
            end

            if isnan(mcmc_use)
                mcmc_use = length(obj.mod.bm_list{1}.samples.s2);
            end

            bassMod = obj.mod.bm_list{1};

            if isnan(prior)
                prior = {};
            end

            p = bassMod.data.p;

            if length(prior) < p
                for i = (length(prior)+1):p
                    tmp.dist = 'uniform';
                    prior{i} = tmp;
                end
            end

            for i = 1:length(prior)
                if ~isfield(prior,'trunc')
                    prior{i}.trunc = [0,1];
                elseif isnan(prior{i}.trunc)
                    prior{i}.trunc = [0,1];
                else
                    prior{i}.trunc = normalizebass(prior{i}.trunc, bassMod.data.bounds(:,i));
                end

                if strcmpi(prior{i}.dist, 'normal') || strcmpi(prior{i}.dist, 'student')
                    prior{i}.mean = normalizebass(prior{i}.mean,bassMod.data.bounds(:,i));
                    prior{i}.sd = prior{i}.sd./(bassMod.data.bounds(2,i)-bassMod.data.bounds(1,i));
                    if strcmpi(prior{i}.dist, 'normal')
                        prior{i}.z = normpdf((prior{i}.trunc(2)-prior{i}.mean)/prior{i}.sd) - normpdf((prior{i}.trunc(1)-prior{i}.mean)/prior{i}.sd);
                    else
                        prior{i}.z = tpdf((prior{i}.trunc(2)-prior{i}.mean)/prior{i}.sd, prior{i}.df) -  tpdf((prior{i}.trunc(1)-prior{i}.mean)/prior{i}.sd, prior{i}.df);
                    end

                    cc = sum(prior{i}.weights.*prior{i}.z);
                    prior{i}.weights = prior{i}.weights./cc;
                end

            end

            pc_mod = obj.mod.bm_list;
            pcs = obj.mod.basis;

            tic
            fprintf('Start\n')

            if int_order > p
                int_order = p;
                warning('int_order > number of inputs, chnage to int_order = number of inputs')
            end

            u_list = arrayfun(@(x) nchoosek(1:p,x), 1:int_order, 'UniformOutput', false);
            ncombs_vec = cellfun(@(x) size(x,1),u_list);
            ncombs = sum(ncombs_vec);
            nxfunc = size(pcs,1);

            n_pc = obj.mod.nbasis;
            w0 = arrayfun(@(x) obj.get_f0(prior,pc_mod,x,mcmc_use), 1:n_pc);

            f0r2 = (pcs*w0').^2;

            tmp = arrayfun(@(x) pc_mod{x}.samples.nbasis(mcmc_use), 1:n_pc);
            max_nbasis = max(tmp);

            C1Basis_array = zeros(n_pc,p,max_nbasis);
            for i = 1:n_pc
                nb = pc_mod{i}.samples.nbasis(mcmc_use);
                mcmc_mod_usei = pc_mod{i}.model_lookup(mcmc_use);
                for j = 1:p
                    for k = 1:nb
                        C1Basis_array(i,j,k) = obj.C1Basis(prior,pc_mod,j,k,i,mcmc_mod_usei);
                    end
                end
            end

            u_list1 = {};
            for i = 1:int_order
                u_list1 = [u_list1; mat2cell(u_list{i},ones(size(u_list{i},1),1),size(u_list{i},2))];
            end

            fprintf('Integrating: %0.2fs\n', toc)

            u_list_temp = [{1:p}; u_list1];

            if ncores > 1
                ints1_temp = cell(length(u_list_temp),1);
                parfor i = 1:length(u_list_temp)
                    ints1_temp{i} = obj.func_hat(prior,u_list_temp{i},pc_mod,pcs,mcmc_use,f0r2,C1Basis_array);
                end
            else
                ints1_temp = arrayfun(@(x) obj.func_hat(prior,x,pc_mod,pcs,mcmc_use,f0r2,C1Basis_array), u_list_temp, 'UniformOutput', false);
            end

            V_tot = ints1_temp{1};
            ints1 = ints1_temp(2:end);

            ints = cell(1,int_order);
            ints{1} = zeros(size(ints1{1},1), size(u_list{1},1));
            for i = 1:size(u_list{1},1)
                ints{1}(:,i) = ints1{i};
            end
            if int_order > 1
                for i = 2:int_order
                    idx = sum(ncombs_vec(1:(i-1)))+(1:size(u_list{i},1));
                    ints{i} = zeros(size(ints1{1},1), length(idx));
                    cnt = 1;
                    for j = idx
                        ints{i}(:,cnt) = ints1{j};
                        cnt = cnt + 1;
                    end
                end
            end

            sob = cell(1,length(u_list));
            sob{1} = ints{1};
            fprintf('Shuffling: %0.2fs\n', toc)

            if length(u_list) > 1
                for i = 2:length(u_list)
                    sob{i} = zeros(nxfunc,size(ints{i},2));
                    for j = 1:size(u_list{i},1)
                        cc = zeros(nxfunc,1);
                        for k = 1:(i-1)
                            ind = arrayfun(@(x) all(ismember(x,u_list{i}(j,:))), u_list{k});
                            cc = cc + (-1)^(i-k)*sum(ints{k}(:,ind),2);
                        end
                        sob{i}(:,j) = ints{i}(:,j)+cc;
                    end
                end
            end

            if isnan(nind)
                nind = ncombs;
            end

            sob_comb_var = cat(2, sob{:});

            vv = mean(sob_comb_var,1);
            [~,ord] = sort(vv,'descend');
            cutoff = vv(ord(nind));
            if nind > length(ord)
                cutoff = min(vv);
            end
            use = sort(find(vv>=cutoff));

            V_other = V_tot - sum(sob_comb_var(:,use),2);
            use = [use ncombs+1];

            sob_comb_var = [sob_comb_var V_other]';
            sob_comb = ((sob_comb_var')./V_tot)';

            sob_comb_var = sob_comb_var(use,:);
            sob_comb = sob_comb(use,:);

            % Calculate "Total Sobol' Index"
            sob_comb_tot = zeros(p, nxfunc);
            idx = 1;
            for i = 1:int_order
                for j = 1:size(u_list{i},1)
                    sob_comb_tot(u_list{i}(j,:), :) = sob_comb_tot(u_list{i}(j,:), :) + sob_comb_var(idx,:);
                    idx = idx + 1;
                end
            end

            names_ind1 = cell(1,ncombs+1);
            cnt = 1;
            for i = 1:length(u_list)
                for j = 1:size(u_list{i},1)
                    tmp = num2str(u_list{i}(j,:));
                    if length(tmp) == 1
                        names_ind1{cnt} = tmp;
                    else
                        tmp = strsplit(tmp);
                        names_ind1{cnt} = strjoin(tmp,'x');
                    end
                    cnt = cnt + 1;
                end
            end
            names_ind1{cnt} = 'other';
            names_ind1 = names_ind1(use);

            fprintf('Finish: %0.2fs\n', toc)

            obj.S = sob_comb;
            obj.S_var = sob_comb_var;
            obj.Var_tot = V_tot;
            obj.T_var = sob_comb_tot;
            obj.names_ind = names_ind1;
            obj.xx = linspace(0,1,nxfunc);

            if ncores > 1
                if isempty(gcp('nocreate'))
                    delete(gcp('nocreate'))
                end
            end

        end

        function out = get_f0(obj, prior, pc_mod, pc, mcmc_use)
            mcmc_mod_use = pc_mod{pc}.model_lookup(mcmc_use);
            out = pc_mod{pc}.samples.beta(mcmc_use,1);
            if (pc_mod{pc}.samples.nbasis(mcmc_use) > 0)
                for m = 1:pc_mod{pc}.samples.nbasis(mcmc_use)
                    out1 = pc_mod{pc}.samples.beta(mcmc_use, 1+m);
                    for l = 1:pc_mod{pc}.data.p
                        out1 = out1.*obj.C1Basis(prior,pc_mod,l,m,pc, mcmc_mod_use);
                    end
                    out = out + out1;
                end
            end
        end

        function out = C1Basis(obj, prior, pc_mod, l, m, pc, mcmc_mod_use)
            int_use_l = find(pc_mod{pc}.samples.vs(mcmc_mod_use,m,:) == l);
            if isempty(int_use_l)
                out = 1;
                return
            end
            s = pc_mod{pc}.samples.signs(mcmc_mod_use,m,int_use_l);
            t = pc_mod{pc}.samples.knots(mcmc_mod_use,m,int_use_l);
            q = 1;

            if s == 0
                out = 0;
                return
            end

            cc = const(s, t);

            if s == 1
                a = max(prior{l}.trunc(1),t);
                b = prior{l}.trunc(2);
                if b < t
                    out = 0;
                    return
                end
                out = obj.intabq1(prior{l},a,b,t,q)/cc;
            else
                a = prior{l}.trunc(1);
                b = min(prior{l}.trunc(2),t);
                if t < a
                    out = 0;
                    return
                end
                out = obj.intabq1(prior{l},a,b,t,q)*(-1)^q/cc;
            end
        end

        function out = C2Basis(obj, prior,pc_mod,l,m1,m2,pc1,pc2,mcmc_mod_use1,mcmc_mod_use2)
            if (l<=pc_mod{pc1}.data.p)
                int_use_l1 = find(pc_mod{pc1}.samples.vs(mcmc_mod_use1,m1,:)==l);
                int_use_l2 = find(pc_mod{pc2}.samples.vs(mcmc_mod_use2,m2,:)==l);

                if isempty(int_use_l1) && isempty(int_use_l2)
                    out = 1;
                    return
                end

                if isempty(int_use_l1)
                    out = obj.C1Basis(prior,pc_mod,l,m2,pc2,mcmc_mod_use2);
                    return
                end

                if isempty(int_use_l2)
                    out = obj.C1Basis(prior,pc_mod,l,m1,pc1,mcmc_mod_use1);
                    return
                end

                q = 1; 
                s1 = pc_mod{pc1}.samples.signs(mcmc_mod_use1,m1,int_use_l1);
                s2 = pc_mod{pc2}.samples.signs(mcmc_mod_use2,m2,int_use_l2);
                t1 = pc_mod{pc1}.samples.knots(mcmc_mod_use1,m1,int_use_l1);
                t2 = pc_mod{pc2}.samples.knots(mcmc_mod_use2,m2,int_use_l2);

                if (t2 < t1)
                    temp = t1;
                    t1 = t2;
                    t2 = temp;
                    temp = s1;
                    s1 = s2;
                    s2 = temp;
                end

                out = obj.C22Basis(prior{l},t1,t2,s1,s2,q);
            else
                % todo categorical
                error('categorical not implemented')
            end
        end

        function out = C22Basis(obj,prior,t1,t2,s1,s2,q)
            cc = const([s1, s2], [t1, t2]);
            if (s1*s2) == 0
                out = 0;
                return
            end

            if (s1 == 1)
                if (s2 == 1)
                    out = obj.intabq2(prior,t2,1,t1,t2,q)./cc;
                    return
                else
                    out = obj.intabq2(prior,t1,t2,t1,t2,q).*(-1)^q/cc;
                    return
                end
            else
                if (s2 == 1)
                    out = 0;
                    return
                else
                    out = obj.intabq2(prior,0,t1,t1,t2,q)/cc;
                    return
                end
            end
        end

        function out = intabq1(obj, prior, a, b, t, q)
            if strcmpi(prior.dist, 'normal')
                if q~=1
                    error('degree other than 1 not supported for normal priors')
                end
                out = 0;
                for k = 1:length(prior.weights)
                    zk = normpdf(b,prior.mean(k),prior.sd(k)) - normpdf(a,prior.mean(k), prior.sd(k));
                    ast = (a-prior.mean(k))/prior.sd(k);
                    bst = (b-prior.mean(k))/prior.sd(k);
                    dnb = normcdf(bst);
                    dna = normcdf(ast);
                    tnorm_mean_zk = prior.mean(k).*zk - prior.sd(k).*(dnb-dna);
                    out = out + prior.weights(k) *(tnorm_mean_zk - t.*zk);
                end
            end
            if strcmpi(prior.dist, 'student')
                if q~=1
                    error('degree other than 1 not supported for student priors')
                end
                out = 0;
                for k = 1:length(prior.weights)
                    int = obj.intx1Student(b,prior.mean(k),prior.sd(k),prior.df(k),t) - intx1Student(a,prior.mean(k),prior.sd(k),prior.df(k),t);
                    out = out + prior.weights(k)*int;
                end
            end
            if strcmpi(prior.dist, 'uniform')
                out = 1/(q+1)*((b-t)^(q+1)-(a-t)^(q+1)) * 1/(prior.trunc(2)-prior.trunc(1));
            end
        end

        function out = intx1Student(obj,x,m,s,v,t)
            temp = (s^2*v)/(m^2 + s^2*v - 2*m*x + x^2);
            out = -((v/(v + (m - x)^2/s^2))^(v/2) * ...
                sqrt(temp) * ...
                sqrt(1/temp) * ...
                (s^2*v* (sqrt(1/temp) - ...
                (1/temp)^(v/2)) + ...
                (t-m)*(-1 + v)*(-m + x) * ...
                (1/temp)^(v/2) * ...
                obj.robust2f1(1/2,(1 + v)/2,3/2,-(m - x)^2/(s^2 *v)) )) / ...
                (s *(-1 + v)* sqrt(v) *beta(v/2, 1/2));
        end

        function out = robust2f1(~,a,b,c,x)
            if abs(x) < 1
                [z,~] = hypergeometric2F1ODE(a,b,c,[0,x]);
                out = z(end);
            else
                [z,~] = hypergeometric2F1ODE(a,c-b,c,0,(1-1/(1-x))/(1-x)^a);
                out = z(end);
            end
        end

        function out = intabq2(obj, prior, a, b, t1, t2, q)
            if strcmpi(prior.dist, 'normal')
                if q~=1
                    error('degree other than 1 not supported for normal priors')
                end
                out = 0;
                for k = 1:length(prior.weights)
                    zk = normpdf(b,prior.mean(k),prior.sd(k)) - normpdf(a,prior.mean(k), prior.sd(k));
                    if (zk < eps)
                        continue
                    end
                    ast = (a-prior.mean(k))/prior.sd(k);
                    bst = (b-prior.mean(k))/prior.sd(k);
                    dnb = normcdf(bst);
                    dna = normcdf(ast);
                    tnorm_mean_zk = prior.mean(k).*zk - prior.sd(k).*(dnb-dna);
                    tnorm_var_zk = zk*prior.sd(k).^2.*(1 + (ast.*dna-bst.*dnb)./zk - ((dna-dnb)./zk).^2) + tnorm_mean_zk.^2./zk;
                    out = out + prior.weights(k) *(tnorm_var_zk - (t1+t2).*tnorm_mean_zk + t1.*t2.*zk);
                    if (out < 0 && abs(out)<1e-12)
                        out = 0;
                    end
                end

            end
            if strcmpi(prior.dist, 'student')
                if q~=1
                    error('degree other than 1 not supported for student priors')
                end
                out = 0;
                for k = 1:length(prior.weights)
                    int = obj.intx2Student(b,prior.mean(k),prior.sd(k),prior.df(k),t1,t2) - obj.intx2Student(a,prior.mean(k),prior.sd(k),prior.df(k),t1,t2);
                    out = out + prior.weights(k)*int;
                end
            end
            if strcmpi(prior.dist, 'uniform')
                out = (sum(obj.pCoef(0:q,q).*(b-t1).^(q-(0:q)).*(b-t2).^(q+1+(0:q))) - sum(obj.pCoef(0:q,q).*(a-t1).^(q-(0:q)).*(a-t2).^(q+1+(0:q)))) .* 1/(prior.trunc(2)-prior.trunc(1));
            end
        end

        function out = intx2Student(~,x,m,s,v,t1,t2)
            temp = (s.^2.*v)/(m.^2 + s.^2.*v - 2*m.*x + x.^2);
            out = ((v/(v + (m - x).^2/s^2)).^(v/2) * ...
                sqrt(temp) * ...
                sqrt(1/temp) * ...
                (-3*(-t1-t2+2*m)*s^2*v* (sqrt(1./temp) - ...
                (1/temp).^(v/2)) + ...
                3*(-t1+m)*(-t2+m)*(-1 + v)*(-m + x) * ...
                (1/temp)^(v/2) * ...
                robust2f1(1/2,(1 + v)/2,3/2,-(m - x)^2/(s^2 *v)) + ...
                (-1+v)*(-m+x)^3*(1/temp)^(v/2) * ...
                robust2f1(3/2,(1 + v)/2,5/2,-(m - x)^2/(s^2 *v)) )) / ...
                (3*s *(-1 + v)* sqrt(v) *beta(v/2, 1/2)); ...

        end

        function out = pCoef(~,i,q)
            out = factorial(q)^2.*(-1).^i./(factorial(q-i).*factorial(q+1+i));
        end

        function out = func_hat(obj,prior,u,pc_mod,pcs,mcmc_use,f0r2,C1Basis_array)
            res = zeros(size(pcs,1),1);
            n_pc = length(pc_mod);
            for i = 1:n_pc
                res = res + pcs(:,i).^2.*obj.Ccross(prior,pc_mod,i,i,u,mcmc_use,C1Basis_array);

                if i < n_pc
                    for j = (i+1):n_pc
                        res = res + 2*pcs(:,i).*pcs(:,j).*obj.Ccross(prior,pc_mod,i,j,u,mcmc_use,C1Basis_array);
                    end
                end
            end

            out = res-f0r2;
        end

        function out = Ccross(obj,prior,pc_mod,i,j,u,mcmc_use,C1Basis_array)
            p = pc_mod{1}.data.p;
            mcmc_mod_usei = pc_mod{i}.model_lookup(mcmc_use);
            mcmc_mod_usej = pc_mod{j}.model_lookup(mcmc_use);

            Mi = pc_mod{i}.samples.nbasis(mcmc_use);
            Mj = pc_mod{j}.samples.nbasis(mcmc_use);

            a0i = pc_mod{i}.samples.beta(mcmc_use,1);
            a0j = pc_mod{j}.samples.beta(mcmc_use,1);
            f0i = obj.get_f0(prior,pc_mod,i,mcmc_use);
            f0j = obj.get_f0(prior,pc_mod,j,mcmc_use);

            out = a0i*a0j + a0i*(f0j-a0j) + a0j*(f0i-a0i);

            if (Mi > 0 && Mj > 0)
                ai = pc_mod{i}.samples.beta(mcmc_use,1+(1:Mi));
                aj = pc_mod{j}.samples.beta(mcmc_use,1+(1:Mj));

                for mi = 1:Mi
                    for mj = 1:Mj
                        temp1 = ai(mi)*aj(mj);
                        temp2 = 1;
                        temp3 = 1;
                        idx = 1:p;
                        if isequal(class(u),'cell')
                            idx2 = u{1};
                        else
                            idx2 = u;
                        end
                        idx(idx2) = [];

                        for l = idx
                            temp2 = temp2.*C1Basis_array(i,l,mi).*C1Basis_array(j,l,mj);
                        end

                        for l = idx2
                            temp3 = temp3.*obj.C2Basis(prior,pc_mod,l,mi,mj,i,j,mcmc_mod_usei,mcmc_mod_usej);
                        end

                        out = out + temp1.*temp2.*temp3;
                    end
                end
            end
        end

        function plot(obj,text,options)
            arguments
                obj
                text = false
                options.labels = [];
                options.col = 'Paired'
                options.time = [];
            end
            col = options.col;
            labels = options.labels;
            time = options.time;
            if isempty(time)
                time = obj.xx;
            end
            if isempty(labels)
                labels1 = obj.names_ind;
            else
                labels1 = obj.names_ind;
                for i = 1:length(labels)
                    labels1{i} = labels{str2double(labels1{i})};
                end
            end
            [map,~,~,~] = brewermap(length(labels1)-1,col);
            rgb = zeros(size(map,1)+1,3);
            rgb(1:size(map,1),:) = map;
            rgb(end,:) = [153,153,153]./255;

            

            [~,ord] = sort(time);
            x_mean = obj.S;
            sens = cumsum(x_mean,1)';
            figure()
            subplot(1,2,1)
            hold all
            idx = find(sum(sens)/size(sens,1)>=.99999,1,'first');
            cnt = 1;
            for i = 1:idx
                x2 = [time(ord), fliplr(time(ord))];
                if i == 1
                    inBetween = [zeros(length(time(ord)),1)', fliplr(sens(ord,i)')];
                else
                    inBetween = [sens(ord,i-1)', fliplr(sens(ord,i)')];
                end
                if mod(cnt,size(rgb,1)+1) == 0
                    cnt = 1;
                end
                fill(x2, inBetween, rgb(cnt,:));
                cnt = cnt + 1;
            end
            ylabel('proportion variance')
            xlabel('x')
            title('Sensitivity')
            ylim([0,1])
            xlim([min(time) max(time)])
            if text
                [~,lab_x] = max(x_mean,[],2);
                cs = zeros(size(sens,2)+1, size(sens,1));
                cs(2:end,:) = cumsum(x_mean,1);
                cs_diff = zeros(size(x_mean,1),size(x_mean,2));
                for i = 1:size(x_mean,2)
                    cs_diff(:,i) = diff(cumsum([0; x_mean(:,1)]));
                end
                tmp = [(1:length(lab_x))' lab_x];
                ind = sub2ind(size(cs),tmp(:,1),tmp(:,2));
                ind1 = sub2ind(size(cs_diff),tmp(:,1),tmp(:,2));
                cs_diff2 = cs_diff./2;
                text(time(lab_x),cs(ind) + cs_diff2(ind1),obj.names_ind)
            end

            subplot(1,2,2)
            x_mean_var = obj.S_var;
            sens_var = cumsum(x_mean_var,1)';
            hold on
            cnt = 1;
            for i = 1:idx
                x2 = [time(ord), fliplr(time(ord))];
                if i == 1
                    inBetween = [zeros(length(time(ord)),1)', fliplr(sens_var(ord,i)')];
                else
                    inBetween = [sens_var(ord,i-1)', fliplr(sens_var(ord,i)')];
                end
                if mod(cnt,size(rgb,1)+1) == 0
                    cnt = 1;
                end
                fill(x2, inBetween, rgb(cnt,:));
                cnt = cnt + 1;
            end
            ylabel('variance')
            xlabel('x')
            title('Variance Decomposition')
            xlim([min(time) max(time)])
            if ~text
                legend(labels1(1:idx),'Location','northwest')
            end

        end
    end
end
