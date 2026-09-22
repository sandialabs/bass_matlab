classdef DmwnchBassTest < BassTestCase
    % Tests for dmwnchBass.m (Walenius' noncentral hypergeometric density,
    % some variables fixed). Deterministic - no RNG involved.

    methods (Test)
        function twoVariableCaseMatchesHandComputation(tc)
            % z_vec = [1;1;2], vars_use = [1 2]:
            %   z_rm = 2, alpha = [0.5; 0.5], j = 2
            %   ss = 1 + 1/(sum(alpha)+1) - sum(1./(alpha+1))
            %      = 1 + 1/2 - (1/1.5 + 1/1.5) = 1.5 - 1.33333... = 1/6
            z_vec = [1; 1; 2];
            vars_use = [1 2];
            out = dmwnchBass(z_vec, vars_use);
            tc.verifyEqual(out, 1/6, 'AbsTol', 1e-12);
        end

        function returnsScalar(tc)
            out = dmwnchBass([1; 2; 3; 4], [1 3]);
            tc.verifySize(out, [1 1]);
            tc.verifyTrue(isfinite(out));
        end

        function probabilityInUnitInterval(tc)
            % Across several fixed-variable subsets the density stays in (0,1].
            z_vec = [2; 3; 1; 4; 2];
            for k = 2:4
                out = dmwnchBass(z_vec, 1:k);
                tc.verifyGreaterThan(out, 0);
                tc.verifyLessThanOrEqual(out, 1 + 1e-12);
            end
        end
    end
end
