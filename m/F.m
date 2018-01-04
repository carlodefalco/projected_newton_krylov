function r = F (x)
  r          = zeros (size (x));
  r(1)       = x(1)^2-1;
  r(2:end-1) = x(1:end-2) - x(2:end-1).^3;
  r(end)     = x(end-1) - x(end);
endfunction 
