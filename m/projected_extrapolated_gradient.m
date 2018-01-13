function [x, err, ll] = projected_extrapolated_gradient (J, F, x0, bounds, opts)
  
  P        = @(x) min (max (x, bounds (:, 1)), bounds (:, 2));

  x        = x0;
  n        = numel (x);

  maxit    = opts.maxit;
  erank    = opts.erank;
  errtol   = opts.errtol;
  maxdamp  = opts.maxdamp;
  lambda0  = opts.lambda0;
  t        = opts.t;
  sigma    = opts.sigma;

  res     = F (x);
  err=normres = norm (res); 
  ll=lambda   = 1;

  for in = 1 : maxit
    
    if (normres <= errtol)
      break
    endif

    err(in) = normres; 
    ll(in)  = lambda;

    semilogy (err)
    drawnow
    
    jac     = J (x);
    X(:, 1) = x;

    for ie = 2 : erank
      d = - jac' * res;

      for m = 0 : maxdamp 
        lambda = lambda0 ^ m;
        xnew = P (x + lambda * d);
        resnew = F (xnew);
        normresnew = norm (resnew);
        if ((1/2) * normresnew^2  <=
            ((1/2) * normres^2 +
             sigma * ((-d)' * (xnew - x)))
            || lambda * (norm (d, inf)) <=
               10 * eps (norm (x, inf)));
          break
        endif
      endfor

      x        = xnew;
      X(:, ie) = x;
      res      = resnew;
      normres  = normresnew;
      jac      = J (x);
      
    endfor

    x       = P (rre (X));
    res     = F (x);
    normres = norm (res);
    
  endfor

  err = err(:);
  ll  = ll(:);
  
endfunction

function s = rre (X)
  
  ## Compute first and second variations
  U = X(:,2:end) - X(:,1:end-1);
  V = U(:,2:end) - U(:,1:end-1);

  ## Eliminate unused u_k column
  U(:,end) = [];
  
  s = X(:,1) - U * (V \ U(:,1));
  
endfunction
