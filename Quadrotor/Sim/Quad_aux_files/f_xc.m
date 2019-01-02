function f_ctrl = f_xc(xc)
%#codegen
g = 9.81;

b_T =  [sin(xc(9)); -cos(xc(9))*sin(xc(8)); cos(xc(9))*cos(xc(8))];

f_ctrl =     [xc(4:6);
              [0;0;g] - xc(7)*b_T;
              zeros(4,1)];
end