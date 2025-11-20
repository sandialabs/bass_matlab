function out = unnormalizebass(z, bounds)
% inverse of normalize

mtmp = repmat(bounds(:,1),1,size(z,1))';
diff = (bounds(:,2) - bounds(:,1));
dtmp = repmat(diff,1,size(z,1))';
out = z .* dtmp + mtmp;
