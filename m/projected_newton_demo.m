clear all
close all

[J, F, x0, bounds] = benchmark_problems ("chenetal");
[x, err, mm, ee, ff, ll] = projected_newton (J, F, x0, bounds);

figure
plot (x)

figure
semilogy (err)
