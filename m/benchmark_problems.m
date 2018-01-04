function [jac, res, x0, bounds] = benchmark_problems (pname)

  switch pname
    case "chenetal"
      jac = @(x) J(x);
      res = @(x) F(x);
      x0          = .9*ones(100,1);
      x0(21:end)  = .5;      
      bounds      = [.5*ones(size (x0)), 2*ones(size (x0))];
      bounds(1,1) = .8;
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

