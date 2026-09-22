classdef BassDataTest < BassTestCase
    % Tests for the BassData.m constructor

    methods (Test)
        function storesDimensionsAndOriginal(tc)
            rng(0, 'twister');
            xx = rand(40, 3);
            y = randn(40, 1);
            bd = BassData(xx, y);

            tc.verifyEqual(bd.n, 40);
            tc.verifyEqual(bd.p, 3);
            tc.verifyEqual(bd.xx_orig, xx);
            tc.verifyEqual(bd.y, y);
        end

        function sumOfSquaresOfResponse(tc)
            y = [1; 2; 3; 4];
            bd = BassData(rand(4, 2), y);
            tc.verifyEqual(bd.ssy, sum(y .^ 2), 'AbsTol', 1e-12);   % 30
        end

        function boundsArePerColumnMinMax(tc)
            xx = [0 10; 5 20; 2 30];
            bd = BassData(xx, randn(3, 1));
            expected = [0 5; 10 30];
            tc.verifyEqual(bd.bounds, expected, 'AbsTol', 1e-12);
        end

        function normalizedInputsInUnitInterval(tc)
            rng(1, 'twister');
            xx = rand(25, 4) * 8 - 3;
            bd = BassData(xx, randn(25, 1));
            tc.verifyEqual(bd.xx, normalizebass(xx, bd.bounds), 'AbsTol', 1e-12);
            tc.verifyGreaterThanOrEqual(bd.xx(:), -1e-12);
            tc.verifyLessThanOrEqual(bd.xx(:), 1 + 1e-12);
        end
    end
end
