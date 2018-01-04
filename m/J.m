function JJ = J (x)
  n   = numel (x);
  d   = [2*x(1); -3*x(2:end-1).^2; -1];
  JJ  = sparse (n,n);
  dm  = [ones(n-1, 1); 0];
  JJ  = spdiags ([dm, d], -1:0, n, n);
endfunction 
