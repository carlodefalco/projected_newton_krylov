

## Copyright (C) 2017 Carlo de Falco

clear all
close all
pkg load bim

KT  = 1e-7;
Lx  = 0;
Rx  = 45;
mu  = 1;
nu  = 2;
g   = 30;
PM  = 30;

hm = @(u) 2 ./  ((1 ./ u(2:end)) + (1 ./ u(1:end-1)));
diffm = @(m,n) hm (mu * g * ((g + 1)/g) .* m .* (n + m) .^ (g-1));
diffn = @(m,n) hm (nu * g * ((g + 1)/g) .* n .* (n + m) .^ (g-1));

G0  = @(p) (200/pi) * atan (4*(PM - p));
G   = @(p) ifelse (G0(p) >= 0, G0(p), 0*p);

pressure = @(m, n) ((g + 1)/g) .* (n + m) .^ g;
dpdm     = @(m, n) (g + 1) .* (n + m) .^ (g - 1);
dpdn     = @(m, n) (g + 1) .* (n + m) .^ (g - 1);

dG0dm    = @(m, n, p) (200/pi) * (-4 * dpdm (m, n)) ./ (1 + (4*(PM - p)) .^2 );;
dGdm     = @(m, n, p) ifelse (G0(p) >= 0, dG0dm (m, n, p), 0*p);

dG0dn    = @(m, n, p) (200/pi) * (-4 * dpdn (m, n)) ./ (1 + (4*(PM - p)) .^2 );;
dGdn     = @(m, n, p) ifelse (G0(p) >= 0, dG0dn (m, n, p), 0*p);

N     = 600;
x     = linspace (Lx, Rx, N+1) .';



dnodes = [];
inodes = [1:N+1];

%m  =  .1 * exp (- (x).^2 / (2 * 5^2));
%n  = 1.1 * exp (- (x - 25).^2 / (2 * 5^2));

m = max (0, .5*sin (2*pi*5*x/45));
n = max (0, .5*cos (2*pi*13*x/45));

mass   = bim1a_reaction (x, ones (N, 1), ones (N+1, 1));

stiffm = @(m, n) bim1a_laplacian (x, diffm (m,n), 1);
stiffn = @(m, n) bim1a_laplacian (x, diffn (m,n), 1);

function [res, jac] = fun (u, mold, nold, inodes, dnodes,
                           stiffm, stiffn, mass,
                           pressure, G, dGdm, dGdn, dt)

  mest = mold;
  nest = nold;
  mest(inodes) = u(1:numel(inodes));
  nest(inodes) = u(numel(inodes)+1:end);

  Sm  = stiffm (mest, nest);
  Sn  = stiffn (mest, nest);
  pp  = pressure (mest, nest);
  GG  = G (pp);  
  Amm = (mass + dt * (Sm - mass .* sparse (diag (GG))));
  Amn = dt * Sm;  
  bm  = mass * mold;
  
  Ann = (mass + dt * Sn);
  Anm = dt * Sn;
  bn  = mass * nold;

  A = [Amm(inodes, inodes), Amn(inodes, inodes);
       Anm(inodes, inodes), Ann(inodes, inodes)];

  b = [(bm(inodes) - Amm(inodes,dnodes) * mold(dnodes, :) -
        Amn(inodes,dnodes) * nold(dnodes, :));
       (bn(inodes) - Anm(inodes,dnodes) * mold(dnodes, :) -
        Ann(inodes,dnodes) * nold(dnodes, :))];
  
  res = A*u-b;

  cfm = mest .* dGdm (mest, nest, pp);
  cfn = mest .* dGdn (mest, nest, pp);
  GGm = mass .* sparse (diag (cfm));
  GGn = mass .* sparse (diag (cfn));
  Amm += (-dt) * GGm;
  Amn += (-dt) * GGn;
  jac = [Amm(inodes, inodes), Amn(inodes, inodes);
         Anm(inodes, inodes), Ann(inodes, inodes)];
  
endfunction

function res = mm (J, h)
  wJ  = J .* h;
  res = (wJ(2:end)+wJ(1:end-1))./(h(2:end)+h(1:end-1));
endfunction

function err = zz (x, u)
  h   = diff (x);
  J   = diff (u) ./ h;
  Jzz = zeros (size (x));
  Jzz(2:end-1) = mm (J, h);
  Jzz(1) = sum (J(1:2) .* h(1:2)) ./ sum (h(1:2));
  Jzz(end) = sum (J(end-1:end) .* h(end-1:end)) ./ sum (h(end-1:end));
  err = zeros (size (J));
  err = ((J - (Jzz(2:end) + Jzz(1:end-1)) / 2) .* h) .^ 2;
endfunction


nt    = 100;
t     = 0;
T     = .05;
tsave = linspace (t, T, nt);
ii    = 2;
tstore= t;
dt    = 1e-5;

msave = mvold = mold = m;
nsave = nvold = nold = n;

for its = 2 : nt
  
  while (t < tsave(its))

    if (t + dt > tsave(its))
      dt = tsave(its) - t;      
    endif
    
    printf ("t = %g, dt = %g\n", t, dt)

    mvold = mold;
    nvold = nold;
    
    mold = m;
    nold = n;

    uguess = [mold;nold];
    if (ii > 3)
      mguess = mvold * (t+dt - t)/(tstore(ii-2) - t) + ...
               mold * (t+dt - tstore(ii-2))/(t - tstore(ii-2));

      nguess = nvold * (t+dt - t)/(tstore(ii-2) - t) + ...
               nold * (t+dt - tstore(ii-2))/(t - tstore(ii-2));

      if (! all ([mguess;nguess] >= 0))
        printf ("negative guess\n")
      endif

      uguess = max (0, [mguess;nguess]);
      
    endif
    
    while (true)

      uold = [mold;nold]; 
      u = fsolve (@(u) ...
                   fun (u, mold, nold, inodes, dnodes,
                        stiffm, stiffn, mass, pressure,
                        G, dGdm, dGdn, dt),
                  uguess, optimset ("Jacobian", "on"));

      all_nzero = all (u >= 0);
      incr = norm ((uguess - u)./(u + 1), inf);
      small_err = incr < 1e-4;

      if (all_nzero && small_err)
        break
      else
        dt *= .5;
        printf (["reducing time step : t = %g, ", ...
                 "dt = %g, all_nzero = %d, small_err = %d\n"],
                t, dt, all_nzero, small_err)
        uguess = [mold;nold];      
        if (ii > 3)
          mguess = mvold * (t+dt - t)/(tstore(ii-2) - t) + ...
                   mold * (t+dt - tstore(ii-2))/(t - tstore(ii-2));

          nguess = nvold * (t+dt - t)/(tstore(ii-2) - t) + ...
                   nold * (t+dt - tstore(ii-2))/(t - tstore(ii-2));

          if (! all ([mguess;nguess] >= 0))
            printf ("negative guess\n")
          endif
          uguess = max (0, [mguess;nguess]);
          
        endif

      endif
      
    endwhile

    
    t += dt;
    tstore(ii) = t;
    if (incr == 0)
      dt *= 2;
    else
      dtfact = min (sqrt(.38) * sqrt (1e-4 / incr), 2);
      dt = min (dtfact * dt, min (diff (tsave)) / 10);
    endif
    
    m = mold;
    m(inodes) = u(1:numel(inodes));

    n = nold;
    n(inodes) = u(numel(inodes)+1:end);

    assert (all (n >= 0))
    assert (all (m >= 0))
    
   
    ii += 1;

  endwhile
  msave (:, end+1) = m;
  nsave (:, end+1) = n;
  figure (1)
  plot (x, m, x, mold, x, n, x, nold)
  title (sprintf ("t = %g", t))
  figure (2)
  plot (x, pressure (m, n))
    
  drawnow
endfor


%%dt m - am div (dm (Dm + Dn)) = G m
