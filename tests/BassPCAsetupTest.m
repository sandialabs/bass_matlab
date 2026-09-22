classdef BassPCAsetupTest < BassTestCase
    % Tests for BassPCAsetup.m (PCA/SVD decomposition of functional response)

    methods (Test)
        function reconstructsCenteredResponse(tc)
            rng(0, 'twister');
            n = 30; q = 8;
            y = randn(n, q);
            setup = BassPCAsetup(y, true, false);

            % basis * newy reconstructs the (transposed) scaled response.
            recon = setup.basis * setup.newy;
            tc.verifyEqual(recon, setup.y_scale', 'AbsTol', 1e-8);
        end

        function eigenvaluesNonNegativeAndDescending(tc)
            rng(1, 'twister');
            y = randn(25, 6);
            setup = BassPCAsetup(y, true, false);
            tc.verifyGreaterThanOrEqual(setup.evals, -1e-10);
            tc.verifyGreaterThanOrEqual(-diff(setup.evals), -1e-10);   % non-increasing
        end

        function centeringSubtractsColumnMeans(tc)
            rng(2, 'twister');
            y = randn(20, 5) + 3;
            setup = BassPCAsetup(y, true, false);
            tc.verifyEqual(setup.y_mean, mean(y, 1), 'AbsTol', 1e-12);
            tc.verifyEqual(mean(setup.y_scale, 1), zeros(1, size(y, 2)), 'AbsTol', 1e-10);
        end

        function noCenteringLeavesMeanAtZero(tc)
            rng(3, 'twister');
            y = randn(15, 4) + 2;
            setup = BassPCAsetup(y, false, false);
            tc.verifyEqual(setup.y_mean, 0);
            tc.verifyEqual(setup.y_scale, y, 'AbsTol', 1e-12);
        end

        function scalingSetsStdAndGuardsZeroVariance(tc)
            rng(4, 'twister');
            y = [randn(20, 3), ones(20, 1)];   % last column has zero variance
            setup = BassPCAsetup(y, true, true);
            tc.verifySize(setup.y_sd, [1, 4]);
            % zero-variance column guarded to 1 (not divided by zero)
            tc.verifyEqual(setup.y_sd(4), 1, 'AbsTol', 1e-12);
            tc.verifyTrue(all(isfinite(setup.y_scale(:))));
        end
    end
end
