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
            bm = bass(obj.x, obj.y(i,:)', obj.opts.nmcmc, obj.opts.nburn, obj.opts.thin, ...
                obj.opts.w1, obj.opts.w2, obj.opts.maxInt, obj.opts.maxBasis, ...
                obj.opts.npart, obj.opts.g1, obj.opts.g2, obj.opts.s2_lower, obj.opts.h1, ...
                obj.opts.h2, obj.opts.a_tau, obj.opts.b_tau, obj.opts.verbose);  % @todo: will need to expand accordingly
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
