%% PROVA FINALE

close all
clear all

n=100;
nn=20;
t = 10^(-4);
sigma = 10^(-4);
x0 =0.9*ones(nn,1);
x0=[x0; 0.5*ones(n-nn,1)]; 
err= norm(F(x0));
verr=err;
tol=10^(-12);
kmax=100;
k=0;
x=x0;
x00=x0;
FLAG=0;
s=zeros(n,1);

nu_k=0.5; 
s0=s;
mmax=20;
vincoli(1,1:2)=[0.8, 2];
vincoli(2:n,1:2 )=[0.5*ones(n-1,1),2*ones(n-1,1)];
vtol=zeros(n,1);
nu_max= 0.9;
sen_1 = zeros(kmax,1);
sen_2= zeros(kmax,1);
lambda0n=0.5;
lambda0g=0.8;   
mm= mmax*ones(kmax, 1);  
P = @(x) min (max (x, vincoli (:, 1)), vincoli (:, 2));

tolleranza_1 = 0.5;
tolleranza_2 = 10^(-1);
count1=0;
count2=0;
count3=0;
epsilon=0.1;
ss=1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load parametri_paper
load XX %sono le soluzioni del paper x del sistema nonlineare ad ogni iterazione
load nu_kv %nu_k ricavate dal paper
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ss=1;
nu_k_vera(end+1)=nu_k_vera(end);
while err>tol && k<22
     k=k+1;
     
     FLAG=direzione_paper(k);
     if FLAG==0 
         
         [d, ~, ~, it, ~, duvec] = gmresx (J(x), -F(x),100 , nu_k_vera(ss)+0.00001, 100, [], [], zeros (100, 1));
        
          ss=ss+1;
          m=mm_paper(k+1);
          lambdap = lambda0n^m;
          lambda(k)=lambdap;
          x0=x;
          x=P(x + lambdap*d);
          sen_1(k)=1;
          mm(k)=m;
         
          
     else
           d= -DOMEGA(F(x),J(x));
           m=mm_paper(k+1);
           lambdap = lambda0g^m;  
           lambda(k)=lambdap;
            x0=x;
            x=P(x + lambdap*d); 
            sen_2(k)=1;
            mm(k)=m;
                 
         end
       err= norm(F(x));
       verr(k+1)=err;
       s=d;
       XXmio(:,end+1)=x;
   
     end

     

sen_1=sen_1(1:k);
sen_2=sen_2(1:k);
mm=mm(1:k);

figure()
ii1=1;
ii2=1;
psen1=[];
psen2=[];
for i=1:k
    if sen_1(i)>0
       psen1(ii1)=i;
       ii1=ii1+1;
       end
    if sen_2(i)>0
       psen2(ii2)=i;
       ii2=ii2+1;
       end
    end


plot(verr, 'k')
ylim([0,max(verr)+1])
hold on 
plot(psen1, verr(psen1),'ob' )
plot(psen2, verr(psen2), 'sr')
plot(verr_paper, '-og')
xlabel('ITERAZIONI')
ylabel('ERRORE')
hold off

%
%end
%end

disp('Il vettore di nu_k per cui funzione è : ')
nu_k_vera+0.00001