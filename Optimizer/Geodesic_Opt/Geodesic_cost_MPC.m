function J = Geodesic_cost_MPC(vars,w,n,K,W_fnc,Phi,Phi_dot)

J = 0;

global GEO_X_MPC;
global GEO_MXDOT_MPC;

for k = 1:K+1 %0 ---> K
    GEO_X_MPC(:,k) = Phi(:,:,k)*vars;
end

w_poly_eval = W_fnc.w_poly_fnc(GEO_X_MPC);

for k = 1:K+1 %0 ---> K
    x_dot_k = Phi_dot(:,:,k)*vars;
    
    W = W_fnc.W_eval(w_poly_eval(:,k));
    M = W\eye(n);
    
    M_xdot = M*x_dot_k;
    J = J + (1/2)*w(k)*(x_dot_k'*M_xdot); 
    
    GEO_MXDOT_MPC(:,k) = M_xdot;
end


return