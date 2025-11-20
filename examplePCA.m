% example using multivariate data
clc; clear

f = @(x) 10. .* sin(pi .* linspace(0,1,50) .* x(:,1)) + 20. .* (x(:,2) - .5).^2 + 10 .* x(:,3) + 5. .* x(:,4);

n = 500;
p = 9;
x = rand(n,p);
xx = rand(1000,p);
y = f(x);

mod = bassPCA(x, y, 5, 99.99, 1);
mod.plot()
pred = mod.predict(xx, [1,100], true);


tmp = squeeze(mean(mod.predict(xx),1))'-f(xx);
disp(var(tmp(:)))

obj = sobolBasis(mod);
obj = obj.decomp(1,NaN,1000,NaN,2);
obj.plot()