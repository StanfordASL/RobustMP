%% Setup MP problem

geodesic_N = 3;
%initializes geodesic_struct for initial invariant set constraint in MP 
setup_geodesic_MPC(n,geodesic_N,W_fnc,dW_fnc,n_W);

%cost matrices
Q = zeros(n); R = eye(m);

%final state constraint
P = 2.5*eye(n);
alpha = 1e-3;

%MP numerics
Tp = 7; %time-horizon
dt = 0.001; %time-resolution
N_mp = 70; %#collocation points
MP_t = (0:dt:Tp)';

%Setup global motion planning problem
[MP_Prob,L_e_mp] = setup_MP(n,m,...
    f,B,df,dB,x_con,u_con,...
    N_mp,Tp,dt,...
    P,alpha,(0.98*d_bar)^2,...
    Q,x_eq,R,u_eq,obs);

%Test MP Solver
if exist('MP_WARM_FLR.mat','file')~=2
    disp('Warning: no warm-start found for MP. Solver may take a while...');
    mp_warm = struct('sol',0);
else
    disp('Loaded MP warm-start');
    load MP_WARM_FLR.mat; %existing stored solution
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
    save('MP_WARM_FLR.mat','mp_warm');
end

%% Setup CCM controller           
[geo_Prob,cntrl_info] = setup_CCM_controller(n,m,geodesic_N,...
                                          W_fnc,dW_fnc,n_W,f,B,...
                                          lambda,MP_state(1,:)',MP_ctrl(1,:)',test_state);