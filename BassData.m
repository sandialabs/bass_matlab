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
        ord     % ord(:, j): row indices sorting xx(:, j)
        rnk     % rnk(i, j): position of row i in ord(:, j)
    end

    methods
        function obj = BassData(xx, y)
            obj.xx_orig = xx;
            obj.y = y;
            obj.ssy = sum(y .* y);
            obj.n = size(xx,1);
            obj.p = size(xx,2);
            obj.bounds = zeros(obj.p, 2);
            for i = 1:obj.p
                obj.bounds(i, 1) = min(xx(:, i));
                obj.bounds(i, 2) = max(xx(:, i));
            end
            obj.xx = normalizebass(obj.xx_orig, obj.bounds);
            [~, obj.ord] = sort(obj.xx, 1);
            obj.rnk = zeros(obj.n, obj.p);
            for j = 1:obj.p
                obj.rnk(obj.ord(:, j), j) = (1:obj.n)';
            end
        end
    end
end
