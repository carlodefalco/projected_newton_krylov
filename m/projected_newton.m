function [x, err, mm, ee, ff, ll] = projected_newton (J, F, x0, bounds)

  maxit    = 2000;
  errtol   = 1e-12;
  maxdamp  = 20;
  
  P        = @(x) min (max (x, bounds (:, 1)), bounds (:, 2));

  x        = x0;
  n        = numel (x);


  lambda0  = .5;
  lambda0G = .8;
  gamma    = .8;
  etamax   = .9;
  eta      =  7.65518617913987e-01;
  eta0     =  7.65518617913987e-01;
  alpha    =  (1 + sqrt (5)) / 2;
  t        =  1.0e-4;
  sigma    =  1.0e-1;
  flag     =  false;

  resnew = F (x);
  normresnew = norm (resnew); 
  
  for in = 1 : maxit
    
    res     = resnew;
    normres = normresnew;
    err(in) = normres; 
    jac     = J (x);
    printf ("%d ", in)
    if (! flag)
      ff(in) = flag;
      if (normres <= errtol)
        ff(in) = flag;
        mm(in) = 0;
        break
      endif
      
      if (in > 1)
        eta = gamma * (err(in) / err(in-1)) ^ alpha;
        maxeta = gamma * ee(in-1) ^ alpha;
        if (maxeta > .1)
          eta = max (eta, maxeta);
        endif
      endif
      %%eta = eta0 * min (1/(in+1)^alpha, normres/(2*err(1)));
      ee(in) = eta;
            
      d = gmres (jac, - res, n, eta, n, [], [], zeros (size (x)));      

      flag = true;
      for m = 0 : maxdamp
        lambda = lambda0 ^ m;
        xnew = P (x + lambda * d);
        resnew = F (xnew); 
        normresnew = norm (resnew);
        printf ("\tPN %g %g %g\n", normresnew,
                (1 - t * lambda * (1-eta)), normres)
        if (normresnew <= ((1 - t * lambda * (1-eta)) * normres))
          flag = false;
          x = xnew;
          break
        endif
      endfor
      ll(in) = lambda;
    endif

    if (flag)
      ff(in) = flag;
      d = - jac' * res;
      flag = true;
      for m = 0 : maxdamp 
        lambda = lambda0G ^ m;
        xnew = P (x + lambda * d);
        resnew = F (xnew);
        normresnew = norm (resnew);
        printf ("\tPG %g %g %g %g\n", 
                ((1/2) * normresnew^2), ...
                ((1/2) * normres^2), ...
                (sigma * (d' * (x - xnew))), ...
                ((1/2) * normres^2 -
                 sigma * (d' * (x - xnew))))
        if ((1/2) * normresnew^2  <=
            ((1/2) * normres^2 -
             sigma * (d' * (x - xnew)))
            || lambda * (norm (d, inf)) <=
               10 * eps (norm (x, inf)));
          flag = false;
          x = xnew;
          break
        endif
      endfor
            
      ee(in) = eta;
      ll(in) = lambda;
    endif

    plot (err)
    drawnow
    mm(in) = m;

  endfor

  err = err(:);
  mm  = mm(:);
  ee  = ee(:);
  ff  = ff(:);
  ll  = ll(:);
  
endfunction
