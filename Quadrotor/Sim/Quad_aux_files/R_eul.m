function R_eul = R_eul(q)
%#codegen

%Body rates to Euler rates
R_eul =  [cos(q(3))/cos(q(2)), -sin(q(3))/cos(q(2)), 0;
          sin(q(3)), cos(q(3)), 0;
         -cos(q(3))*tan(q(2)), sin(q(3))*tan(q(2)), 1];
          
end