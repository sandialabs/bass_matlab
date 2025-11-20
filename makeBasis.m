function temp2 = makeBasis(signs, vs, knots, xdata)
% Make basis funciton using continuous variables

cc = const(signs, knots);

a = bsxfun(@times, signs(:)', bsxfun(@minus, xdata(:, vs), knots(:)'));
temp1 = max(0,a);
if length(signs) == 1
    temp2 = temp1/cc;
else
    temp2 = prod(temp1, 2)./cc;
end