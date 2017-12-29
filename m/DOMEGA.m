function ddO=DOMEGA(f,J)
    n=max(size(f));
    ddO= zeros (n,1);
    
    
    ddO(1,1)= f(1)*(J(1,1)) + f(2);
    for i=2:n-1
        ddO(i,1)= f(i)*(J(i,i)) + f(i+1); 
            
        end    
    ddO(n,1)= -f(n);
    
    end