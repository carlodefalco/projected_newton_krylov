function [jac, res, x0, bounds, opts] = benchmark_problems (pname)

  switch pname
    case "chenetal_paper"
      jac = @(x) J(x);
      res = @(x) F(x);
      x0          = .9*ones(100,1);
      x0(21:end)  = .5;      
      bounds      = [.5*ones(size (x0)), 2*ones(size (x0))];
      bounds(1,1) = .8;
      opts.maxit    = 400;
      opts.errtol   = 1e-12;
      opts.maxdamp  = 20;
      opts.lambda0  = .5;
      opts.lambda0G = .8;
      opts.gamma    =  1;
      opts.etamax   = .8;
      opts.eta0     =  7.65518617913987e-01;
      opts.alpha    =  2;
      opts.t        =  1.0e-4;
      opts.sigma    = -2.575;  %% THIS VALUE IS NOT ALLOWED ACCORDING
                               %% TO THE PAPER BUT IS REQUIRED TO
                               %% ACHIEVE CONVERGENCE
    case "chenetal"
      jac = @(x) J(x);
      res = @(x) F(x);
      x0          = .9*ones(100,1);
      x0(21:end)  = .5;      
      bounds      = [.5*ones(size (x0)), 2*ones(size (x0))];
      bounds(1,1) = .8;
      opts.maxit    = 400;
      opts.errtol   = 1e-12;
      opts.maxdamp  = 20;
      opts.lambda0  = .5;
      opts.lambda0G = .8;
      opts.gamma    =  0.9;
      opts.etamax   = .9;
      opts.eta0     =  7.65518617913987e-01;
      opts.alpha    =  2;
      opts.t        =  1.0e-4;
      opts.sigma    =  1.0e-4;  
  
    case "diffreactmonotone"
      pkg load bim msh fpl
      n = 102;
      x = linspace (0, 1, n+2);
      A = bim1a_laplacian (x, 1, 1);
      M = bim1a_reaction  (x, 1, 1);
      jac = @(x) Jdr (x, n, A, M);
      res = @(x) Fdr (x, n, A, M);
      x0 = zeros (n, 1);
      bounds      = [zeros(size (x0)), 2*ones(size (x0))];
      opts.maxit    = 100;
      opts.errtol   = 1e-12;
      opts.maxdamp  = 20;
      opts.lambda0  = .5;
      opts.lambda0G = .8;
      opts.gamma    =  1;
      opts.etamax   = .8;
      opts.eta0     = .8;
      opts.alpha    =  2;
      opts.t        =  1.0e-4;
      opts.sigma    =  1.0e-4;
      
    case "scalarlocmin"
      x = [0, 1, 2, 3, 4];
      y = [-1, 1, 1e-2, 1, -1];
      p = polyfit (x, y, 4);
      dp = polyder (p);
      jac = @(xx) polyval (dp, xx);
      res = @(xx) polyval (p, xx);
      x0 = 2.1;
      bounds  = [0, 3];
      opts.maxit    = 20;
      opts.errtol   = 1e-12;
      opts.maxdamp  = 20;
      opts.lambda0  = .5;
      opts.lambda0G = .8;
      opts.gamma    =  1;
      opts.etamax   = .8;
      opts.eta0     = .8;
      opts.alpha    =  2;
      opts.t        =  1.0e-4;
      opts.sigma    =  1.0e-4;
    otherwise
      error ("unknown problem name")
  endswitch
  
endfunction

function JJ = J (x)
  n   = numel (x);
  d   = [2*x(1); -3*x(2:end-1).^2; -1];
  JJ  = sparse (n,n);
  dm  = [ones(n-1, 1); 0];
  JJ  = spdiags ([dm, d], -1:0, n, n);
endfunction 

function r = F (x)
  n          = numel (x);
  r          = zeros (n, 1);
  r(1)       = x(1)^2-1;
  r(2:end-1) = x(1:end-2) - x(2:end-1).^3;
  r(end)     = x(end-1) - x(end);
endfunction 

function JJ = Jdr (x, n, A, M)
  JJ = A(2:end-1, 2:end-1) + M(2:end-1, 2:end-1) .* diag (exp (x) + exp (-x));
endfunction 

function r = Fdr (x, n, A, M)
  r = A(2:end-1, 2:end-1) * x + M(2:end-1, 2:end-1) * (exp (x) - exp (-x)) - M(2:end-1, 2:end-1) * ones (size (x));
endfunction 

