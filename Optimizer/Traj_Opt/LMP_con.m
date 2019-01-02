function c = LMP_con(xu,Prob,n,m,N,P,D,f,B,Tp,MP_t,MP_x)

global Geo_vel_init geodesic_MPC;

%actual initial state
x_init = Prob.user.x_act;

no = Prob.user.obs.n_obs;

c = zeros(n*(N+1)+2+no*(N+1),1);

%% Dynamics constraints

c(1:n*(N+1)) = (2/Tp)*D*xu(1:n*(N+1)) -...
                dyn_all(f,B,xu(1:n*(N+1)),xu(n*(N+1)+1:end-1),n,m,N);

%% Initial RPI constraint

[~, X_dot,J_opt,~,geo_result,~] = compute_geodesic(geodesic_MPC.geo_Prob,...
    n,geodesic_MPC.geodesic_N,xu(1:n),x_init,geodesic_MPC.T_e,geodesic_MPC.T_dot_e,geodesic_MPC.warm,geodesic_MPC.solver);
Geo_vel_init = X_dot(:,1);
geodesic_MPC.warm.sol = 1; geodesic_MPC.warm.result = geo_result;

c(n*(N+1)+1) = J_opt;

%% Terminal constraint

T_star = xu(end);
x_star = interp1(MP_t,MP_x,T_star);

c(n*(N+1)+2) = (xu(1+N*n:(N+1)*n)-x_star')'*P*(xu(1+N*n:(N+1)*n)-x_star');

%% Obstacle constraints

for i = 1:no
    o_pos = Prob.user.obs.pos(:,i);
    for k = 1:N+1
        Mo = Prob.user.obs.M_obs(:,:,i,k);
        x_k = xu(1+(k-1)*n:2+(k-1)*n);
        c(n*(N+1)+2+(i-1)*(N+1)+k,1) = (o_pos-x_k)'*Mo*(o_pos-x_k);
    end
end

end

