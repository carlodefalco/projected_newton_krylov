function y = F(x)
    x = x(:);
    y = zeros(size(x));
%    y(1)=x(1)^2-x(2)-2;
%    y(2)=x(1)-x(2);
    
    y(1)= x(1)^2-1;
    n=numel(x);
    for i=2:n-1
        y(i)=x(i-1)-x(i)^3;
        end
    y(n)=x(n-1)-x(n);
    end