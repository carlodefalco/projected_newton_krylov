clc
clear all
nn=20;
n=100;
clear R
  
%Carico:
load parametri_paper
% verr_paper: la norma del residuo delle iterazioni del paper ;
% mm_paper : le m del paper;
% direzione_paper: i valori che assume la variabile FLAG nel paper (0 Newton-K, 1 direzione del gradiente) 


function y = F(x)
    x = x(:);
    y = zeros(size(x));
    y(1)= x(1)^2-1;
    n=numel(x);
    for i=2:n-1
        y(i)=x(i-1)-x(i)^3;
        end
    y(n)=x(n-1)-x(n);
    end


function M = J(x)
      n=numel(x);
      M=sparse(n, n);
      M(1,1)=2*x(1);
      
      for i=2:n-1
          M(i,i-1)=1;
          M(i,i)=-3*x(i)^2;
          end
      M(n,n-1)=1;
      M(n,n)=-1;
    end
    

function ddO=DOMEGA(f,J)
    ddO=J'*f;
    end    

    
x0 =0.9*ones(nn,1);
x0=[x0; 0.5*ones(n-nn,1)]; 
k=2; % !! parte da seconda iterazione del paper!!

vincoli(1,1:2)=[0.8, 2];
vincoli(2:n,1:2 )=[0.5*ones(n-1,1),2*ones(n-1,1)];

P = @(x) min (max (x, vincoli (:, 1)), vincoli (:, 2));


while k<24
      printf('ITERAZIONE %d\n', k)
      
      if (k~=6 && k~=10 && k~=14 )
          lambdap= 0.5^(mm_paper(k));

          
          r0=F(x0);
          J0=J(x0);
          
          [d, ~, ~, it, ~, duvec] = gmresx (J0, -r0, 100 , 0, 100, [], [], zeros (100, 1));
          clear R;
          jj=-1;
          for j=1:columns(duvec)
            
             R(j)=norm(F(P(x0+lambdap*duvec(:,j))));
             
             if k < 17
               soglia=0.00002;
             end
             
             if k >= 17
               soglia=0.001;
             end
             
             
             if abs(R(j)-verr_paper(k))/R(j)<soglia
                 jj=j;
                 break;
             endif   
             
          end
          
%%% QUI HO SELEZIONATO IO I VALORI DELL'ERRORE CHE SONO PIÙ VICINI A QUELLI DEL PAPER
%%% PERCHÈ SONO TROPPO PICCOLI PER LA SOGLIA STABILITA

           if k ==20
            jj=6;
            end
           if k ==22
            jj=18;
            end
          if k ==23
            jj=2;
            end
          disp('valore dell errore nel paper: ')
          verr_paper(k)
          disp('posizione dell errore in R : ')
          jj
          disp('i primi valori di R : ')
          R
           jjj(end+1)=jj;
          if (jj==-1)
            disp('NON È STATO TROVATO UN ERRORE UGUALE!!')
            break;
          endif
          
          err1(end+1)=norm(J(x0)*duvec(:,jj)+F(x0));
          err2(end+1)=norm(J(x0)*duvec(:,jj-1)+F(x0));
          err3(end+1)=norm(J(x0)*duvec(:,jj+1)+F(x0));
          nu_k_vera(end+1)=err1(end)/verr_paper(k-1);
          nu_k1(end+1)=err1(end)/verr_paper(k-1);
          nu_k2(end+1)=err2(end)/verr_paper(k-1);
          nu_k3(end+1)=err3(end)/verr_paper(k-1);
          
          x=P(x0 + lambdap*duvec(:,jj));
          x0=x;
          % XX contiene le soluzioni x di ogni iterazione k
          XX(:,end+1)=x;
          k=k+1;
          DD(:,end+1)=duvec(:,jj);
      
      
      else
      
          dd0=DOMEGA(F(x0), J(x0));
          lambdap =0.8^mm_paper(k);% sul paper c'è 0.64, ma in realtà sunziona solo con 0.8
         
          printf('DIREZIONE DEL GRADIENTE CON LAMBDA = %d\n',lambdap);
          x=P(x0 - lambdap*dd0);
          verr_paper(k)=norm(F(x));
          norm(F(x))
          XX(:,end+1)=x;
          x0=x; 
          k=k+1; 
          

      end
endwhile    
%nu_k_vera=(nu_k1+nu_k2)/2;
%nu_k_vera=(nu_k2+nu_k3)/2;

save XX XX
save nu_kv nu_k_vera
save jjj jjj
save DD DD