function Qf = getQf(XtX, Xty)
% Get the quadratic form y'X solve(X'X) X'y,
% as well as least squares beta and cholesky of X'X

try
    R = chol(XtX);
catch
    Qf.fullrank = false;
    return
end

dr = diag(R);
if length(dr) > 1
    if max(dr(2:end))/min(dr) > 1e3
        Qf.fullrank = false;
        return
    end
end

tmp1 = (R')\Xty;
bhat = R\tmp1;
qf = bhat' * Xty;
Qf.R = R;
Qf.bhat = bhat;
Qf.qf = qf;
Qf.fullrank = true;
