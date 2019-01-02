%% Setup (global) MP problem

geodesic_N = 3; %Chebyshev order of geodesic 
%initializes geodesic_struct for initial geodesic constraint when
%generating motion plans 
setup_geodesic_MPC(n,geodesic_N,W_fnc,dW_fnc,n_W);

%cost matrices
Q = 0*diag([1;1;0;0;0;0]); R = eye(m);

%final state constraint: \|x-x_eq\|_P^2 \leq \alpha
P = 2.5*eye(n);
alpha = 1e-3;

%MP numerics
Tp = 19; %horizon
N_mp = 120; %#collocation points
dt = 0.001; %time resolution
MP_t = (0:dt:Tp)';

%Setup global motion planning problem
[MP_Prob,L_e_mp] = setup_MP(n,m,...
    f,B,df,dB,x_con,u_con,...
    N_mp,Tp,dt,...
    P,alpha,(0.98*d_bar)^2,...
    Q,x_eq,R,u_eq,obs);

%Test MP Solver
if exist('MP_WARM_PVTOL.mat','file')~=2
    disp('Warning: no warm-start found for MP. Solver may not converge with full obstacle set');
    mp_warm = struct('sol',0);
else
    disp('Loaded MP warm-start');
    load MP_WARM_PVTOL.mat; %existing stored solution
end

% MP Solve
tic
[MP_state,MP_ctrl,converged_MP,mp_warm] = compute_MP(MP_Prob,...
    test_state,x_con,u_con,x_eq,u_eq,...
    n,m,N_mp,L_e_mp,mp_warm);
toc
disp('MP:'); disp(converged_MP);

if (converged_MP >= 1 && converged_MP <= 3)
    %store successful solution
    mp_warm.sol = 1.0;
    save('MP_WARM_PVTOL.mat','mp_warm');
end

%% setup (local) MP problem - uses time-varying tubes

if (do_mpc)
    
    T_mpc = 3; %MPC horizon
    delta = 1; %MPC resolve time
    N_mpc = 18; %#collocation points
    
    mu = 4; %weight for re-join time
    
    %adjust obstacles for MPC
    obs_mpc = obs;
    M_obs_mpc = zeros(2,2,obs.n_obs,N_mpc+1);
    for i = 1:obs.n_obs
        M_obs_mpc(:,:,i,:) = repmat(M_obs(:,:,i),1,1,N_mpc+1);
    end
    obs_mpc.M_obs = M_obs_mpc;
    
    %Setup unscaled ellipsoids (for online dynamically scaled obstacles)
    obs_mpc.time_var = time_var; %to enable online re-sized tubes
    [U_u,S_u,V_u] = svd(M_ccm_pos_unscaled);
    obs_mpc.U = U_u; obs_mpc.V = V_u; obs_mpc.S = S_u;
    obs_mpc.r = obs_rad + len;
    
    %Setup MPC problem
    [MPC_Prob,L_e_mpc,MPC_st] = setup_LMP(MP_state,MP_ctrl,MP_t,...
        n,m,...
        f,B,df,dB,x_con,u_con,...
        N_mpc,T_mpc,delta,dt,...
        P,alpha,d_bar^2,...
        Q,x_eq,R,u_eq,mu,obs_mpc);
    
    %Test MPC solver
    tic
    mpc_warm = struct('sol',0,'Tp',T_mpc,'s_t',MPC_st,'d_bar',d_bar);
    [MPC_state,~,~,converged_MPC,mpc_warm,MPC_Prob,obs_active] = compute_LMP(MPC_Prob,...
        test_state,x_con,u_con,MP_t,MP_state,MP_ctrl,...
        n,m,N_mpc,L_e_mpc,mpc_warm,dt,0,delta,(d_bar)^2,lambda,obs_mpc);
    toc
    disp('MPC:');disp(converged_MPC);
    
    MPC_Prob.CHECK = 1;
    mpc_warm.sol = 1;
    
end

%% Setup CCM controller           
[geo_Prob,cntrl_info] = setup_CCM_controller(n,m,geodesic_N,...
                                          W_fnc,dW_fnc,n_W,f,B,...
                                          lambda,MP_state(1,:)',MP_ctrl(1,:)',test_state);