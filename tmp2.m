function p = makeParabola(a,b,c)
   p = @parabola;

   function y = parabola(x)
   y = a*x.^2 + b*x + c;
   end

end