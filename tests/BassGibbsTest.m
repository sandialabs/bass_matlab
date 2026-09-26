classdef BassGibbsTest < BassTestCase
    % Tests for the Gibbs-block options (intercept, center_basis, s2_fixed)
    % and the bassGibbs stepping class.

    methods (TestMethodSetup)
        function seed(~)
            rng(2026, 'twister');
        end
    end

    methods (Test)
        function explicitDefaultsRun(tc)
            [X, y] = BassGibbsTest.makeData(100);
            mod = bass(X, y, 'nmcmc', 200, 'nburn', 100, 'verbose', false, ...
                'intercept', true, 'center_basis', false);
            tc.verifyTrue(all(isfinite(mod.predict(X)), 'all'));
        end

        function centeredNoInterceptIsMeanZero(tc)
            % every posterior draw averages to zero over the training inputs
            [X, y] = BassGibbsTest.makeData(100);
            y = y - mean(y);
            mod = bass(X, y, 'nmcmc', 400, 'nburn', 300, 'verbose', false, ...
                'intercept', false, 'center_basis', true);
            P = mod.predict(X);
            tc.verifyLessThan(max(abs(mean(P, 2))), 1e-10);
            tc.verifyGreaterThan(max(mod.samples.nbasis), 0);
            tc.verifyLessThan(mean((mean(P, 1)' - y).^2), 0.5 * var(y));
        end

        function fixedS2StaysFixed(tc)
            [X, y] = BassGibbsTest.makeData(80);
            mod = bass(X, y, 'nmcmc', 100, 'nburn', 50, 'verbose', false, 's2_fixed', 0.7);
            tc.verifyEqual(mod.samples.s2, 0.7 * ones(50, 1));
        end

        function predictFollowsMcmcUseOrder(tc)
            [X, y] = BassGibbsTest.makeData(80);
            mod = bass(X, y, 'nmcmc', 300, 'nburn', 200, 'verbose', false);
            P = mod.predict(X);
            idx = [50 3 20 3];
            tc.verifyEqual(mod.predict(X, 'mcmc_use', idx), P(idx, :), 'AbsTol', 1e-12);
        end

        function gibbsFittedMatchesPredict(tc)
            [X, y] = BassGibbsTest.makeData(80);
            for cfg = {{true, false}, {false, true}}
                g = bassGibbs(X, y - mean(y), 10, 'g1', 2, 'g2', 0.1, ...
                    'intercept', cfg{1}{1}, 'center_basis', cfg{1}{2});
                mu = zeros(10, numel(y));
                for k = 1:10
                    g.setResponse(y - mean(y) + 0.1 * randn(size(y)));
                    g.step(5);
                    mu(k, :) = g.fitted()';
                    g.save();
                end
                tc.verifyEqual(g.model().predict(X), mu, 'AbsTol', 1e-8);
            end
        end

        function gibbsRecoversLatentRegression(tc)
            % z_i ~ N(f(x_i), s2), y_i ~ N(z_i, 0.1^2): the BASS block tracks f
            [X, ~, f] = BassGibbsTest.makeData(150);
            ftrue = f(X) - mean(f(X));
            y = ftrue + 0.3 * randn(size(ftrue)) + 0.1 * randn(size(ftrue));
            g = bassGibbs(X, y, 50, 'intercept', false, 'center_basis', true, ...
                'g1', 2, 'g2', 0.05);
            z = y;
            for it = 1:400
                g.setResponse(z);
                g.step(3);
                mu = g.fitted(); s2 = g.s2();
                v = 1 / (1 / s2 + 1 / 0.01);
                z = v * (mu / s2 + y / 0.01) + sqrt(v) * randn(size(y));
                if it > 350, g.save(); end
            end
            muhat = mean(g.model().predict(X), 1)';
            tc.verifyLessThan(mean((muhat - ftrue).^2), 0.25 * var(ftrue));
        end
    end

    methods (Static)
        function [X, y, f] = makeData(n)
            f = @(x) 10 * sin(pi * x(:,1) .* x(:,2)) + 20 * (x(:,3) - .5).^2;
            X = rand(n, 3);
            y = f(X) + randn(n, 1);
        end
    end
end
