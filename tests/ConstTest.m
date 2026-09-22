classdef ConstTest < BassTestCase
    % Tests for const.m (max value of a BASS basis function on [0,1] inputs)

    methods (Test)
        function singleFactorPositiveSign(tc)
            % (1+1)/2 - 1*0.5 = 0.5
            tc.verifyEqual(const(1, 0.5), 0.5, 'AbsTol', 1e-12);
        end

        function multipleFactorsProduct(tc)
            % prod((1+1)/2 - 1*0) over two factors = 1*1 = 1
            tc.verifyEqual(const([1; 1], [0; 0]), 1, 'AbsTol', 1e-12);
        end

        function zeroGuardReturnsOne(tc)
            % (-1+1)/2 - (-1)*0 = 0 -> guarded to 1
            tc.verifyEqual(const(-1, 0), 1, 'AbsTol', 1e-12);
        end

        function mixedSignProduct(tc)
            % factor1: (1+1)/2 - 1*0.25 = 0.75
            % factor2: (-1+1)/2 - (-1)*0.75 = 0.75
            % product = 0.5625
            tc.verifyEqual(const([1; -1], [0.25; 0.75]), 0.5625, 'AbsTol', 1e-12);
        end

        function acceptsRowOrColumnVectors(tc)
            % const uses (:) internally, so orientation must not matter
            tc.verifyEqual(const([1 -1], [0.25 0.75]), const([1; -1], [0.25; 0.75]), ...
                'AbsTol', 1e-12);
        end
    end
end
