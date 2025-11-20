classdef BassPrior
    % Structure to store prior

    properties
        maxInt
        maxBasis
        npart
        g1
        g2
        s2_lower
        h1
        h2
        a_tau
        b_tau
        w1
        w2
    end

    methods
        function obj = BassPrior(maxInt, maxBasis, npart, g1, g2, s2_lower, h1, h2, a_tau, b_tau, w1, w2)
            obj.maxInt = maxInt;
            obj.maxBasis = maxBasis;
            obj.npart = npart;
            obj.g1 = g1;
            obj.g2 = g2;
            obj.s2_lower = s2_lower;
            obj.h1 = h1;
            obj.h2 = h2;
            obj.a_tau = a_tau;
            obj.b_tau = b_tau;
            obj.w1 = w1;
            obj.w2 = w2;
        end
    end
end
