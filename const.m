function cc = const(signs, knots)
% Get max value of BASS basis function, assuming 0-1 range of inputs

cc = prod((signs(:)+1)/2-signs(:).*knots(:));

if cc == 0
    cc = 1;
end
