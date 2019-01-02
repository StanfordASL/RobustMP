function xi = phi(xc)
%#codegen

g = 9.81;
b_T =  [sin(xc(9));-cos(xc(9))*sin(xc(8));cos(xc(9))*cos(xc(8))];

xi =  [xc(1:6);
      [0;0;g]-xc(7)*b_T;
       xc(10)];
end