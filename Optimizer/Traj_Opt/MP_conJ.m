function conJ = MP_conJ(xu,Prob,n,m,N,P,D,df,dB,Tp,x_eq,obs)
%Dynamics, and terminal

no = obs.n_obs;

global Geo_vel_init geodesic_MPC B_full;

conJ = zeros(n*(N+1)+2+no*(N+1),(n+m)*(N+1));

%% Dynamics constraints

conJ(1:n*(N+1),1:n*(N+1)) = (2/Tp)*D - ...
            df_all(df,xu(1:n*(N+1)),n,N) - ...
            dB_all(dB,xu(1:n*(N+1)),...
                      xu(n*(N+1)+1:end),n,m,N);

conJ(1:n*(N+1),n*(N+1)+1:end) = -B_full;

%% Initial RPI constraint

w_poly = geodesic_MPC.W.w_poly_fnc(xu(1:n));
M = (geodesic_MPC.W.W_eval(w_poly))\eye(n);
conJ(n*(N+1)+1,1:n) = -2*Geo_vel_init'*M;


%% Terminal constraint
conJ(n*(N+1)+2,n*N+1:n*(N+1)) = 2*(P*(xu(n*N+1:n*(N+1))-x_eq))';

%% Obstacle constraints

for i = 1:no
    o_pos = obs.pos(:,i);
    Mo = obs.M_obs(:,:,i);
    for k = 1:N+1
        x_k = xu(1+(k-1)*n:2+(k-1)*n);
        conJ(n*(N+1)+2+(i-1)*(N+1)+k,1+(k-1)*n:2+(k-1)*n) = -2*(o_pos-x_k)'*Mo';
    end
end

end

