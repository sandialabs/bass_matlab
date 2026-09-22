classdef MakeBasisTest < BassTestCase
    % Tests for makeBasis.m (tensor-product hinge basis functions)

    methods (Test)
        function singleVariableHingeMatchesFormula(tc)
            xdata = [0; 0.25; 0.5; 0.75; 1];
            signs = 1;
            vs = 1;
            knots = 0.5;
            out = makeBasis(signs, vs, knots, xdata);

            cc = const(signs, knots);
            expected = max(0, signs * (xdata(:, vs) - knots)) / cc;
            tc.verifyEqual(out, expected, 'AbsTol', 1e-12);
        end

        function hingeIsZeroBelowKnotAndNonNegative(tc)
            xdata = linspace(0, 1, 11)';
            out = makeBasis(1, 1, 0.5, xdata);
            % values at inputs strictly below the knot are zero
            tc.verifyEqual(out(xdata < 0.5), zeros(sum(xdata < 0.5), 1), 'AbsTol', 1e-12);
            tc.verifyGreaterThanOrEqual(out, 0);
        end

        function returnsColumnVectorMatchingRowCount(tc)
            rng(0, 'twister');
            xdata = rand(20, 3);
            out = makeBasis([1; -1], [1 3], [0.4; 0.6], xdata);
            tc.verifySize(out, [20, 1]);
        end

        function twoVariableTensorProductMatchesHandComputation(tc)
            % One row so we can hand-verify the product form.
            xdata = [0.8, 0.2, 0.3];   % 1x3
            signs = [1; -1];
            vs = [1, 3];
            knots = [0.5; 0.5];
            out = makeBasis(signs, vs, knots, xdata);

            cc = const(signs, knots);
            f1 = max(0, 1 * (0.8 - 0.5));     % 0.3
            f2 = max(0, -1 * (0.3 - 0.5));    % 0.2
            expected = (f1 * f2) / cc;
            tc.verifyEqual(out, expected, 'AbsTol', 1e-12);
        end
    end
end
