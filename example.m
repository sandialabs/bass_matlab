% example using 1-D data
clc; clear

f = @(x) 10 * sin(2*pi * x(:,1) .* x(:,2)) + 20. * (x(:,3) - .5).^2 + ...
    10 * x(:,4) + 5. * x(:,5);
n = 500;
p = 10;
x = rand(n,p);
xx = rand(1000,p);
y = f(x) + randn(n,1);

mod = bass(x, y);
pred = mod.predict(xx, [1,100], true);

mod.plot()

disp(var(mean(mod.predict(xx),1)'-f(xx)))
