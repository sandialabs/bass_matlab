classdef GetQfTest < BassTestCase
    % Tests for getQf.m (quadratic form, Cholesky, least-squares beta)

    methods (Test)
        function fullRankReturnsCorrectQuadraticForm(tc)
            rng(0, 'twister');
            X = randn(50, 4);
            XtX = X' * X;
            y = randn(50, 1);
            Xty = X' * y;

            Qf = getQf(XtX, Xty);

            tc.verifyTrue(Qf.fullrank);
            tc.verifyEqual(Qf.R' * Qf.R, XtX, 'AbsTol', 1e-8, 'RelTol', 1e-10);
            tc.verifyEqual(Qf.bhat, XtX \ Xty, 'AbsTol', 1e-8, 'RelTol', 1e-8);
            tc.verifyEqual(Qf.qf, Xty' * (XtX \ Xty), 'AbsTol', 1e-8, 'RelTol', 1e-8);
        end

        function rankDeficientReturnsNotFullRank(tc)
            % Duplicate column -> XtX is singular, chol fails.
            X = [1 1; 2 2; 3 3];
            XtX = X' * X;
            Xty = X' * [1; 2; 3];
            Qf = getQf(XtX, Xty);
            tc.verifyFalse(Qf.fullrank);
        end

        function illConditionedDiagonalRejected(tc)
            % Diagonal Cholesky factor ratio > 1e3 triggers the guard.
            XtX = diag([1, 1e8]);
            Xty = [1; 1];
            Qf = getQf(XtX, Xty);
            tc.verifyFalse(Qf.fullrank);
        end

        function scalarFullRankCase(tc)
            XtX = 4;
            Xty = 8;
            Qf = getQf(XtX, Xty);
            tc.verifyTrue(Qf.fullrank);
            tc.verifyEqual(Qf.bhat, 2, 'AbsTol', 1e-12);
            tc.verifyEqual(Qf.qf, 16, 'AbsTol', 1e-12);
        end
    end
end
