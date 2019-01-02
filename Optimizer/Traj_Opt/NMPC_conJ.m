function conJ = NMPC_conJ(xu,Prob,n,m,N,P,D,f,df,B,dB,Tp,x_eq)

global Geo_vel_init geodesic_MPC B_full;
global NMPC_CON_J;

%% Dynamics constraints

NMPC_CON_J(1:n*(N+1),1:n*(N+1)) = (2/Tp)*D - ...
            df_all(df,xu(1:n*(N+1)),n,N) - ...
            dB_all(dB,xu(1:n*(N+1)),...
                                xu(n*(N+1)+1:end),n,m,N);

NMPC_CON_J(1:n*(N+1),n*(N+1)+1:end-1) = -B_full;

%% Initial RPI constraint

w_poly = geodesic_MPC.W.w_poly_fnc(xu(1:n));
M = (geodesic_MPC.W.W_eval(w_poly))\eye(n);
NMPC_CON_J(n*(N+1)+1,1:n) = -2*Geo_vel_init'*M;

%% Terminal constraint

NMPC_CON_J(n*(N+1)+2,1+N*n:(N+1)*n) = 2*(xu(1+N*n:(N+1)*n)-x_eq)'*P;

conJ = NMPC_CON_J;

end

