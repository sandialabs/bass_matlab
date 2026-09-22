classdef BassSmokeTest < BassTestCase
    % Seeded end-to-end smoke test of bass() and BassModel.predict().
    % Uses short MCMC and small data - checks plumbing, shapes, and finiteness,
    % NOT statistical accuracy. Exercises BassState.update (birth/death/change),
    % writeState, makeBasisMatrix, genCandBasis, genBasisChange, logProbChangeMod.

    properties
        x
        y
        xx
        opts
    end

    methods (TestMethodSetup)
        function seedAndBuildData(tc)
            rng(0, 'twister');
            f = @(x) 10 * sin(2*pi * x(:,1) .* x(:,2)) + 20 * (x(:,3) - 0.5).^2 ...
                + 10 * x(:,4) + 5 * x(:,5);
            n = 80; p = 5;
            tc.x = rand(n, p);
            tc.y = f(tc.x) + randn(n, 1);
            tc.xx = rand(15, p);
            tc.opts = {'nmcmc', 300, 'nburn', 200, 'thin', 1, 'verbose', false};
        end
    end

    methods (Test)
        function returnsBassModelWithExpectedStoreCount(tc)
            mod = bass(tc.x, tc.y, tc.opts{:});
            tc.verifyClass(mod, 'BassModel');
            tc.verifyEqual(mod.nstore, floor((300 - 200) / 1));
        end

        function predictHasCorrectShapeAndIsFinite(tc)
            mod = bass(tc.x, tc.y, tc.opts{:});
            pred = mod.predict(tc.xx);
            tc.verifySize(pred, [mod.nstore, size(tc.xx, 1)]);
            tc.verifyTrue(all(isfinite(pred(:))));
        end

        function predictSubsetAndNuggetShape(tc)
            mod = bass(tc.x, tc.y, tc.opts{:});
            pred = mod.predict(tc.xx, 'mcmc_use', [1 5], 'nugget', true);
            tc.verifySize(pred, [2, size(tc.xx, 1)]);
            tc.verifyTrue(all(isfinite(pred(:))));
        end

        function thinningReducesStoreCount(tc)
            mod = bass(tc.x, tc.y, 'nmcmc', 300, 'nburn', 200, 'thin', 2, 'verbose', false);
            tc.verifyEqual(mod.nstore, floor((300 - 200) / 2));
            tc.verifyEqual(numel(mod.samples.s2), mod.nstore);
        end
    end
end
