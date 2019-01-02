function c = NMPC_con(xu,Prob,n,m,N,P,D,f,B,Tp,x_eq)

global Geo_vel_init geodesic_MPC;

%actual initial state
x_init = Prob.user.x_act;

c = zeros(n*(N+1)+2,1);

%% Dynamics constraints

c(1:n*(N+1)) = (2/Tp)*D*xu(1:n*(N+1)) -...
                dyn_all(f,B,xu(1:n*(N+1)),xu(n*(N+1)+1:end),n,m,N);

%% Initial RPI constraint

[~, X_dot,J_opt,~,geo_result,~] = compute_geodesic(geodesic_MPC.geo_Prob,...
    n,geodesic_MPC.geodesic_N,xu(1:n),x_init,geodesic_MPC.T_e,geodesic_MPC.T_dot_e,geodesic_MPC.warm,geodesic_MPC.solver);
Geo_vel_init = X_dot(:,1);
geodesic_MPC.warm.sol = 1; geodesic_MPC.warm.result = geo_result;

c(n*(N+1)+1) = J_opt;

%% Terminal constraint

c(n*(N+1)+2) = (xu(1+N*n:(N+1)*n)-x_eq)'*P*(xu(1+N*n:(N+1)*n)-x_eq);

end

