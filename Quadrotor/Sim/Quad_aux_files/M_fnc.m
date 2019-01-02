function M = M_fnc(xc,M_xi)
%#codegen

b_T =  [sin(xc(9)); -cos(xc(9))*sin(xc(8)); cos(xc(9))*cos(xc(8))];
db_T_q =    [0, cos(xc(9));
            -cos(xc(8))*cos(xc(9)), sin(xc(8))*sin(xc(9));
            -sin(xc(8))*cos(xc(9)),-cos(xc(8))*sin(xc(9))];
        
%Jacobian of transform from xc -> xic        
Phi_jac =  blkdiag(eye(3),eye(3),-[b_T, db_T_q*xc(7)],1);

%M(xc) = Phi_jac(xc)'*M_xi*Phi_jac(xc)
%W(xc) = Psi_jac(xic)'*W_xi*Psi_jac(xic)

M =  M_xi*Phi_jac; %leave out first part


end