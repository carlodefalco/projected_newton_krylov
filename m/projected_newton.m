function [x, err, mm, ee] = projected_newton (J, F, x0, bounds)

  P        = @(x) min (max (x, bounds (:, 1)), bounds (:, 2));

  x        = x0;
  n        = numel (x);


  lambda0  = .5;
  gamma    =  1;
  eta      = .5;
  alpha    = 2;
  t        = 1e-4;


  for in = 1 : 200

    res = F(x);
    
    err(in) = norm (res);
    if (err(in) <= 1e-12)
      break
    endif
    jac = J(x);
    
    if (in > 1)
      eta = min (gamma * (err(in) / err(in-1))^alpha, .3);
      maxeta = gamma*ee(in-1)^alpha;
      if (maxeta < .3 && maxeta > .1)
        eta = max (eta, maxeta);
      endif
    endif
    ee(in) = eta;
    
    %%d = jac \ ((- res));
    tol = eta;
    [d, flag, relres] = gmres (jac, - res, n, tol, n, [], [], zeros (size (x)));
    

    for m = 0:36
      lambda = lambda0 ^ m;
      xnew = P(x+lambda*d);
      if (norm (F(xnew)) <= ((1 - t*lambda * (1-eta)) * err(in)))
        break
      endif 
    endfor

    mm(in) = m;
    x = xnew;
    
    %plot (1:n, x, 1:n, ones (n, 1))
    %axis ([1, numel(x), 0, 2])
    %drawnow

  endfor

endfunction
