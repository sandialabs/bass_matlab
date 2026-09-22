classdef BassPCASmokeTest < BassTestCase
    % Seeded end-to-end smoke test of bassPCA() and BassBasis.predict() for
    % functional/multivariate response. Checks plumbing, shapes, and finiteness,
    % NOT statistical accuracy. Covers bassPCA -> BassPCAsetup -> BassBasis ->
    % per-PC bass -> functional predict.

    properties
        x
        y
        xx
        q
        npc
    end

    methods (TestMethodSetup)
        function seedAndBuildData(tc)
            rng(0, 'twister');
            tc.q = 20;
            f = @(x) 10 * sin(pi * linspace(0, 1, tc.q) .* x(:,1)) ...
                + 20 * (x(:,2) - 0.5).^2 + 10 * x(:,3) + 5 * x(:,4);
            n = 80; p = 4;
            tc.x = rand(n, p);
            tc.y = f(tc.x);                 % n x q
            tc.xx = rand(12, p);
            tc.npc = 3;
        end
    end

    methods (Test)
        function returnsBassBasis(tc)
            mod = bassPCA(tc.x, tc.y, tc.npc, 99.99, 1, true, false, ...
                'nmcmc', 300, 'nburn', 200, 'thin', 1, 'verbose', false);
            tc.verifyClass(mod, 'BassBasis');
            tc.verifyEqual(mod.nbasis, tc.npc);
            tc.verifyEqual(numel(mod.bm_list), tc.npc);
        end

        function functionalPredictHasThreeDimShape(tc)
            mod = bassPCA(tc.x, tc.y, tc.npc, 99.99, 1, true, false, ...
                'nmcmc', 300, 'nburn', 200, 'thin', 1, 'verbose', false);
            nstore = mod.bm_list{1}.nstore;
            pred = mod.predict(tc.xx);
            % [npoints x q x nstore]
            tc.verifySize(pred, [size(tc.xx, 1), tc.q, nstore]);
            tc.verifyTrue(all(isfinite(pred(:))));
        end

        function predictSubsetOfMcmcSamples(tc)
            mod = bassPCA(tc.x, tc.y, tc.npc, 99.99, 1, true, false, ...
                'nmcmc', 300, 'nburn', 200, 'thin', 1, 'verbose', false);
            pred = mod.predict(tc.xx, [1 2 3]);
            tc.verifySize(pred, [size(tc.xx, 1), tc.q, 3]);
            tc.verifyTrue(all(isfinite(pred(:))));
        end
    end
end
