clear all
close all

[J, F, x0, bounds] = benchmark_problems ("diffreactmonotone");
[x, err, mm, ee] = projected_newton (J, F, x0, bounds);

figure
plot (x)

figure
semilogy (err)
