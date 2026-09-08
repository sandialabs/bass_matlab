classdef PoolBass
    % class for parallel BASS

    properties
        x
        y
        opts
    end

    methods
        function obj = PoolBass(x, y, opts)
            obj.x = x;
            obj.y = y;
            obj.opts = opts;
            obj.opts.verbose = false;
        end

        function bm = rowbass(obj, i)
            bm = bass(obj.x, obj.y(i,:)', 'nmcmc', obj.opts.nmcmc, 'nburn', obj.opts.nburn, ...
                'thin', obj.opts.thin, 'w1', obj.opts.w1, 'w2', obj.opts.w2, ...
                'maxInt', obj.opts.maxInt, 'maxBasis', obj.opts.maxBasis, ...
                'npart', obj.opts.npart, 'g1', obj.opts.g1, 'g2', obj.opts.g2, ...
                's2_lower', obj.opts.s2_lower, 'h1', obj.opts.h1, 'h2', obj.opts.h2, ...
                'a_tau', obj.opts.a_tau, 'b_tau', obj.opts.b_tau, 'verbose', obj.opts.verbose);
        end

        function out = fit(obj, ncores, nrow_y)
            if isempty(gcp('nocreate'))
                parpool(ncores);
            end
            out = cell(1,nrow_y);
            bar = ProgressBar(nrow_y, ...
                'IsParallel', true, ...
                'WorkerDirectory', pwd, ...
                'Title', 'Running MCMC Chains' ...
                );
            bar.setup([], [], []);
            parfor i = 1:nrow_y
                out{i} = obj.rowbass(i);
                updateParallel([], pwd);
            end
            bar.release();
        end
    end
end
