%% Prova per testare se le eta_k ricavate da "reverse_end.m" funzionano

close all
clear all

n     = 100;
nn    = 20;
t     = 1.0e-4;
sigma = 1.0e-4;
x0    = 0.9 * ones (nn, 1);
x0    = [x0; (0.5 * ones (n - nn, 1))]; 
err   = norm (F (x0));
verr  = err;
tol   = 1.0e-12;
kmax  = 100;
k     = 0;
x     = x0;
x00   = x0;
FLAG  = 0;
s     = zeros (n, 1);

eta_k  = 0.5; 
s0    = s;
mmax  = 20;

vincoli(1,1:2)    = [0.8, 2];
vincoli(2:n,1:2 ) = [(0.5 * ones (n-1, 1)), (2 * ones (n-1, 1))];
vtol     = zeros (n, 1);
eta_max   = 0.9;
sen_1    = zeros (kmax, 1);
sen_2    = zeros (kmax, 1);
lambda0n = 0.5;
lambda0g = 0.8;   
mm = mmax * ones (kmax, 1);  
ss      = 1;

P  = @(x) min (max (x, vincoli(:, 1)), vincoli(:, 2));

%% carico alcuni parametri utili
load parametri_paper
load XX            % sono le soluzioni x del paper del sistema nonlineare ad
                   % ogni iterazione
load eta_kv        % eta_k ricavate dal paper
eta_k_definitive =  eta_k_vera + 0.00001;

%% aggiusto le dimensioni
eta_k_definitive = [eta_k_definitive(1:4), 0, eta_k_definitive(5:8), ...
                   eta_k_definitive(8), eta_k_definitive(9:11), ...
                   eta_k_definitive(11), eta_k_definitive(12:end)];

                   
FLAG = 0;
eta_k_vera(end+1) = eta_k_vera(end);

while ((err > tol) && (k++ < 25))
  
  if (FLAG == 0)
    
    eta_k = eta_k_definitive(ss);
    [d, ~, ~, it, ~, duvec] = gmresx (J(x), -F(x), 100,
                                      eta_k, 100, [], [],
                                      zeros (100, 1));
    ss = ss + 1;
    m  = 0;
    while (m < (mmax+1))

      lambdap = lambda0n ^ m;
      
      if (norm (F(P (x + lambdap * d))) <=
          ((1- lambdap * t * (1 - eta_k)) * norm (F (x))))
        
        x0 = x;
        x = P (x + lambdap * d);
        sen_1(k) = 1;
        FLAG = 0;
        mm(k) = m;
        m = mmax+1;

      else
        
        m = m+1;
        FLAG = 1;
        
      endif
      
      
    endwhile
    lambda(k) = lambdap;
    
  else
    d = -DOMEGA (F(x), J(x));
    lambdap = 0.8;  %qui non sono riuscita ad aggiungere il ciclo while con la 
                    %disuguaglianza (6) perchè, come si vede dai risultati stampati, 
                    % questa non è MAI verificata per lambdap=0.8.
    lambdac = 0.8;
    printf ("iterazione = %d\n", k)

    printf ("lhs della disuguaglianza (6) = %d\n",
            OMEGA(F (P (x + lambdac*d))) - OMEGA (F (x)))

    printf ("rhs della disuguaglianza (6) = %d\n",
            sigma * DOMEGA (F (x), J (x))' *
            (P (x + lambdac * d) - x))
    
    x0 = x;
    x = P (x + lambdap * d); 
    sen_2(k) = 1;
    FLAG = 0;
    mm(k) = m;
    m = mmax + 1;
    lambda(k) = lambdap;
    FLAG = 0;
  endif
  err = norm (F (x));
  
  verr(k+1) = err;
  s = d;
  if (FLAG ~= 1)
    XXmio(:,end+1) = x;
  endif
endwhile


sen_1 = sen_1(1:k);
sen_2 = sen_2(1:k);
mm    = mm(1:k);

figure()
ii1=1;
ii2=1;
psen1=[];
psen2=[];
for i=1:k
  if sen_1(i)>0
    psen1(ii1)=i;
    ii1=ii1+1;
  endif
  if sen_2(i)>0
    psen2(ii2)=i;
    ii2=ii2+1;
  endif
endfor

verr_paper2=[verr_paper(1:4), verr_paper(4), verr_paper(5:9), ...
             verr_paper(9), verr_paper(10:12), verr_paper(13), ...
             verr_paper(13:end)];
plot(verr, 'k')
ylim([0,max(verr)+1])
hold on 
plot(verr_paper2, '-og')

plot(psen1, verr(psen1),'ob' )
plot(psen2, verr(psen2), 'sr')
xlabel('ITERAZIONI')
ylabel('ERRORE')
hold off

rm XX                 
