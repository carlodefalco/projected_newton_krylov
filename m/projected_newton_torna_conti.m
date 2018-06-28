%% è uguale al file projected_newton.m solo che qui 
% ho  imposto lambda = 0.8 ogni volta che 
% si usa la direzione del gradiente. 

function [x, err, mm, ee, ff, ll] = projected_newton_torna_conti (J, F, x0, bounds, opts)
  
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

  resnew = F (x);
  normresnew = norm (resnew); 

  eta      = eta0;
  m        = 0;
  lambda   = lambda0^m;
  flag     = false;
  % carico il vettore di eta_k calcolate dal revers_eng.m
  load eta_kv
  eta_k_definitive = eta_k_vera +0.00001;
  eta_k_definitive = [eta_k_definitive(1:4), 0, eta_k_definitive(5:8), ...
       eta_k_definitive(8), eta_k_definitive(9:11), ...
       eta_k_definitive(11), eta_k_definitive(12:end)];


  ff(1)= 0;
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
    mm(in)  = m;
    
    if (in > 1)
      eta = gamma * (err(in) / err(in-1)) ^ alpha;
      eta = min (eta, etamax);
      %DISABLE "PRACTICAL SAFEGUARD"
      maxeta = gamma * ee(in) ^ alpha;
      if (maxeta > .1)
        eta = max (eta, maxeta);
      endif
       
      
    endif

    
    d = gmres (jac, - res, n, eta, n, [], [], zeros (size (x)));      

    flag = false;
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
x = xnew;
    if (! flag)

      printf ("%10.10s\n", "PN")
      ff (end + 1) = flag;
      
    else

      printf ("%10.10s\n","PG")
      ff (end + 1) = flag;
      
      d = - jac' * res;
      flag = true;
      m = 1;
      for m = 0 : 20  % impongo lambda = 0.8!!! 
        lambda =  lambda0G ^ m;
        xnew = P (x + lambda * d);
        resnew = F (xnew);
        normresnew = norm (resnew);
        if ((1/2) * normresnew^2  <=   %
            ((1/2) * normres^2 +    %
             sigma * ((-d)' * (xnew - x))))  %
%            %|| lambda * (norm (d, inf)) <=
%              % 10 * eps (norm (x, inf)));
          flag = false;
          x = xnew;
          break  %
        endif  %
      endfor
            
    endif

%    plot (err)
%    drawnow

  endfor

  err = err(:);
  mm  = mm(:);
  ee  = ee(:);
  ff  = ff(:);
  ll  = ll(:);
 
 figure ();
 hold on;
 for i = 1:size(err)(1)
    if (ff(i)==0)  
        plot (i, err(i), 'ob');
    elseif    
        plot (i, err(i), 'or');
    endif  
  endfor
  xlabel ('iterations');
  ylabel ('||F(x)||');
  plot (err, '--k');
 endfunction
