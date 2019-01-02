function u_aux = compute_opt_aux(X,X_dot,E,...
                            W_fun,f_fun,B_fun,u_nom,lambda)

%X: n x K_e matrix: [x_1....x_{Ke}] geodesic points
%X_dot: n x K_e matrix: [x_dot_1...x_dot_{Ke}] geodesic curve velocities

%% min norm control

K_e = size(X,2);

w_poly = W_fun.w_poly_fnc(X);
W_start = W_fun.W_eval(w_poly(:,1));
W_end = W_fun.W_eval(w_poly(:,K_e));

A = 2*X_dot(:,K_e)'*(W_end\B_fun(X(:,K_e)));

% l_b = -Inf;
u_b = -2*lambda*E + 2*X_dot(:,1)'*(W_start\(f_fun(X(:,1)) + B_fun(X(:,1))*u_nom)) -...
                    2*X_dot(:,K_e)'*(W_end\(f_fun(X(:,K_e)) + B_fun(X(:,K_e))*u_nom));

a = -u_b;
b = A';
if (norm(b)==0) || (a <= 0) 
    u_aux = zeros(length(b),1);
else
    u_aux = -(a/(b'*b))*b;
end

end



