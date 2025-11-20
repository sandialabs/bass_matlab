classdef BassData
    % Structure to store data

    properties
        xx_orig
        y
        ssy
        n
        p
        bounds
        xx
    end

    methods
        function obj = BassData(xx, y)
            obj.xx_orig = xx;
            obj.y = y;
            obj.ssy = sum(y .* y);
            obj.n = length(xx);
            obj.p = size(xx,2);
            obj.bounds = zeros(obj.p, 2);
            for i = 1:obj.p
                obj.bounds(i, 1) = min(xx(:, i));
                obj.bounds(i, 2) = max(xx(:, i));
            end
            obj.xx = normalizebass(obj.xx_orig, obj.bounds);
        end
    end
end
