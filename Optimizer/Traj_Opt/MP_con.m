function c = MP_con(xu,Prob,n,m,N,P,D,f,B,Tp,x_eq,obs)
%Dynamics, and terminal

%actual initial state
x_init = Prob.user.x_act;

no = obs.n_obs;

global Geo_vel_init geodesic_MPC;

c = zeros(n*(N+1)+2+no*(N+1),1);

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
c(n*(N+1)+2) = (xu(n*N+1:n*(N+1))-x_eq)'*P*(xu(n*N+1:n*(N+1))-x_eq);

%% Obstacle constraints

for i = 1:no
    o_pos = obs.pos(:,i);
    Mo = obs.M_obs(:,:,i);
    for k = 1:N+1
        x_k = xu(1+(k-1)*n:2+(k-1)*n);
        c(n*(N+1)+2+(i-1)*(N+1)+k,1) = (o_pos-x_k)'*Mo*(o_pos-x_k);
    end
end

end

