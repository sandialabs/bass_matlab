function out = normalizebass(x, bounds)
% normalize to 0-1 scale
mtmp = repmat(bounds(:,1),1,size(x,1))';
diff = (bounds(:,2) - bounds(:,1));
dtmp = repmat(diff,1,size(x,1))';

out = (x-mtmp)./dtmp;
