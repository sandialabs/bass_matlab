classdef NormalizeBassTest < BassTestCase
    % Tests for normalizebass.m and unnormalizebass.m

    methods (Test)
        function normalizesColumnwiseToUnitInterval(tc)
            x = [0 10; 5 20; 10 30];          % 3x2
            bounds = [0 10; 10 30];           % per-column [min max]
            out = normalizebass(x, bounds);
            expected = [0 0; 0.5 0.5; 1 1];
            tc.verifyEqual(out, expected, 'AbsTol', 1e-12);
        end

        function unitBoundsLeaveDataUnchanged(tc)
            x = [0.1 0.9; 0.4 0.6];
            bounds = [0 1; 0 1];
            out = normalizebass(x, bounds);
            tc.verifyEqual(out, x, 'AbsTol', 1e-12);
        end

        function roundTripRecoversOriginal(tc)
            rng(0, 'twister');
            x = randn(15, 4) * 3 + 7;
            bounds = [min(x)' max(x)'];
            z = normalizebass(x, bounds);
            back = unnormalizebass(z, bounds);
            tc.verifyEqual(back, x, 'AbsTol', 1e-10);
        end

        function normalizedValuesWithinZeroOne(tc)
            rng(1, 'twister');
            x = rand(30, 3) * 100 - 50;
            bounds = [min(x)' max(x)'];
            z = normalizebass(x, bounds);
            tc.verifyGreaterThanOrEqual(z(:), -1e-12);
            tc.verifyLessThanOrEqual(z(:), 1 + 1e-12);
        end
    end
end
