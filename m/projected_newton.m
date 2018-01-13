function [x, err, mm, ee, ff, ll] = projected_newton (J, F, x0, bounds, opts)
  
  P        = @(x) min (max (x, bounds (:, 1)), bounds (:, 2));

  x        = x0;
  n        = numel (x);

  maxit    = opts.maxit;
  errtol   = opts.errtol;
  maxdamp  = opts.maxdamp;
  lambda0  = opts.lambda0;
  lambda0G = opts.lambda0G;
  gamma    = opts.gamma;
  etamax   = opts.etamax;
  eta0     = opts.eta0;
  alpha    = opts.alpha;
  t        = opts.t;
  sigma    = opts.sigma;

  resnew     = F (x);
  normresnew = norm (resnew); 
  err        = normresnew;
  
  ee=eta      = eta0;
  mm=m        = 0;
  ll=lambda   = lambda0^m;
  ff=flag     = false;


  printf ("%10.10s | \t%10.10s | \t%10.10s | \t%s | \t%10.10s\n",
          "  Iterates", "    lambda",
          "     ||F||", "       eta",
          "      flag");

  printf ("__________________________________________________________________________\n");

  for in = 1 : maxit
    
    res     = resnew;
    normres = normresnew;
    jac     = J (x);

    if (normres <= errtol)
      break
    endif

    printf ("%10.5d |\t%10.5g |\t%10.5g |\t%10.5g |\t",
            in, lambda, normres, eta);

    ee(in)  = eta;
    err(in) = normres; 
    ll(in)  = lambda;
    ff(in)  = flag;
    mm(in)  = m;

    
    if (in > 1)
      eta = gamma * (err(in) / err(in-1)) ^ alpha;
      eta = min (eta, etamax);
      %% DISABLE "PRACTICAL SAFEGUARD"
      %%maxeta = gamma * ee(in) ^ alpha;
      %%if (maxeta > .1)
      %%  eta = max (eta, maxeta);
      %%endif
    endif

    
    d = gmres (jac, - res, n, eta, n, [], [], zeros (size (x)));      

    flag = true;
    for m = 0 : maxdamp
      lambda = lambda0 ^ m;
      xnew = P (x + lambda * d);
      resnew = F (xnew); 
      normresnew = norm (resnew);
      if (normresnew <= ((1 - t * lambda * (1-eta)) * normres))
        flag = false;
        x = xnew;
        break
      endif
    endfor

    if (! flag)

      printf ("%10.10s\n", "PN")
      
    else

      printf ("%10.10s\n","PG")
      
      d = - jac' * res;
      flag = true;
      for m = 0 : maxdamp 
        lambda = lambda0G ^ m;
        xnew = P (x + lambda * d);
        resnew = F (xnew);
        normresnew = norm (resnew);
        if ((1/2) * normresnew^2  <=
            ((1/2) * normres^2 +
             sigma * ((-d)' * (xnew - x)))
            || lambda * (norm (d, inf)) <=
               10 * eps (norm (x, inf)));
          flag = false;
          x = xnew;
          break
        endif
      endfor
            
    endif

    semilogy (err)
    drawnow
    
  endfor

  err = err(:);
  mm  = mm(:);
  ee  = ee(:);
  ff  = ff(:);
  ll  = ll(:);
  
endfunction
