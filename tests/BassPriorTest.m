classdef BassPriorTest < BassTestCase
    % Tests for the BassPrior.m constructor. Guards against regressions in the
    % positional argument order (maxInt, maxBasis, npart, g1, g2, s2_lower,
    % h1, h2, a_tau, b_tau, w1, w2).

    methods (Test)
        function storesAllFieldsInOrder(tc)
            bp = BassPrior(3, 1000, 20, 0, 0, 0, 10, 10, 0.5, 250, 5, 5);

            tc.verifyEqual(bp.maxInt, 3);
            tc.verifyEqual(bp.maxBasis, 1000);
            tc.verifyEqual(bp.npart, 20);
            tc.verifyEqual(bp.g1, 0);
            tc.verifyEqual(bp.g2, 0);
            tc.verifyEqual(bp.s2_lower, 0);
            tc.verifyEqual(bp.h1, 10);
            tc.verifyEqual(bp.h2, 10);
            tc.verifyEqual(bp.a_tau, 0.5);
            tc.verifyEqual(bp.b_tau, 250);
            tc.verifyEqual(bp.w1, 5);
            tc.verifyEqual(bp.w2, 5);
        end

        function distinctValuesAreNotSwapped(tc)
            % Use all-distinct values so a mis-ordering would be caught.
            bp = BassPrior(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12);
            fields = {'maxInt','maxBasis','npart','g1','g2','s2_lower', ...
                'h1','h2','a_tau','b_tau','w1','w2'};
            for i = 1:numel(fields)
                tc.verifyEqual(bp.(fields{i}), i, ...
                    sprintf('field %s should equal %d', fields{i}, i));
            end
        end
    end
end
