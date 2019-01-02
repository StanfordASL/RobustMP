function conJ = LMP_conJ(xu,Prob,n,m,N,P,D,f,df,B,dB,Tp,MP_t,MP_x,MP_u)

obs = Prob.user.obs;
no = obs.n_obs;

global Geo_vel_init geodesic_MPC B_full;
global LMP_CON_J;

%% Dynamics constraints

LMP_CON_J(1:n*(N+1),1:n*(N+1)) = (2/Tp)*D - ...
            df_all(df,xu(1:n*(N+1)),n,N) - ...
            dB_all(dB,xu(1:n*(N+1)),...
                                xu(n*(N+1)+1:end-1),n,m,N);

LMP_CON_J(1:n*(N+1),n*(N+1)+1:end-1) = -B_full;

%% Initial RPI constraint

w_poly = geodesic_MPC.W.w_poly_fnc(xu(1:n));
M = (geodesic_MPC.W.W_eval(w_poly))\eye(n);

LMP_CON_J(n*(N+1)+1,1:n) = -2*Geo_vel_init'*M;

%% Terminal constraint

T_star = xu(end);
x_star = interp1(MP_t,MP_x,T_star);
u_star = interp1(MP_t,MP_u,T_star);

LMP_CON_J(n*(N+1)+2,1+N*n:(N+1)*n) = 2*(xu(1+N*n:(N+1)*n)-x_star')'*P;
LMP_CON_J(n*(N+1)+2,end) = -2*(xu(1+N*n:(N+1)*n)-x_star')'*P*(f(x_star')+B(x_star')*u_star');

%% Obstacle constraints

for i = 1:no
    o_pos = Prob.user.obs.pos(:,i);
    for k = 1:N+1
        Mo = Prob.user.obs.M_obs(:,:,i,k);
        x_k = xu(1+(k-1)*n:2+(k-1)*n);
        LMP_CON_J(n*(N+1)+2+(i-1)*(N+1)+k,1+(k-1)*n:2+(k-1)*n) = -2*(o_pos-x_k)'*Mo';
    end
end

conJ = LMP_CON_J;

end

