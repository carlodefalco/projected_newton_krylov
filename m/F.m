function y = F(x)
    x = x(:);
    y = zeros(size(x));
    y(1) = x(1)^2-1;
    n = numel(x);
    for i=2:n-1
      y(i) = x(i-1)-x(i)^3;
    endfor
    y(n) = x(n-1)-x(n);
endfunction