function m = lomed (x)
%x is the input vector, n is the size of x
x = sort(x);
n=size(x,2);
if rem(n, 2) == 1
    %if n odd
   m = x((n+1)/2);
else
    %if n is even
   m = x(n/2);
end











