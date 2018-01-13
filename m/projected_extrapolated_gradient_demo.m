clear all
close all

jac         = @(x) J(x);
res         = @(x) F(x);
x0          = .9*ones(100,1);
x0(21:end)  = .5;      
bounds      = [.5*ones(size (x0)), 2*ones(size (x0))];
bounds(1,1) = .8;

opts.maxit    = 400;
opts.erank    = 5;
opts.errtol   = 1e-12;
opts.maxdamp  = 20;
opts.lambda0  = .5;
opts.t        =  1.0e-4;
opts.sigma    =  1.0e-4;  


[x, err, ll] = projected_extrapolated_gradient (jac, res, x0, bounds, opts);

figure
plot (x)

figure
semilogy (err)
