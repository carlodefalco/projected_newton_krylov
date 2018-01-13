clear all
close all

benchmarks = {'chenetal', 'chenetal_paper', 'diffreactmonotone', ...
              'scalarlocmin'};
choice = menu ('select benchmark', benchmarks);
[J, F, x0, bounds, opts] = benchmark_problems (benchmarks{choice});
[x, err, mm, ee, ff, ll] = projected_newton_torna_conti(J, F, x0, bounds, opts);

figure
plot (x)

figure
semilogy (err)
