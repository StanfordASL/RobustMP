%% Setup Geodesic numerics

geodesic_N = 1;

setup_geodesic_MPC(n,geodesic_N,W_fnc,dW_fnc,n_W); %initializes geodesic_MPC struct 

%% Setup MPC numerics

%lookahead
T_mpc = 1.5;

%re-solve time
delta = 0.1;

%time resolution of solution
dt = 0.001;

%collocation points
N_mpc = 24;

% Setup MPC problem 
[MPC_Prob,L_e_mpc,MPC_st] = setup_NMPC(n,m,...
    f,B,df,dB, state_constr,ctrl_constr,...
    N_mpc,T_mpc,dt,...
    P,alpha,(d_bar)^2,...
    Q,x_eq,R,u_eq);

%load MPC_WARM_TMPC.mat;
 mpc_warm = struct('Tp',T_mpc,'shift',0,'sol',0,'solve_t',0,...
                   's_t',MPC_st,'state',[],'ctrl',[],'result',[]);

%% Test MPC Solve
      
tic
[MP_state,MP_ctrl,converged_MPC,mpc_warm,MPC_Prob] = compute_NMPC(MPC_Prob,...
    test_state,state_constr,ctrl_constr,x_eq,u_eq,...
    n,m,N_mpc,L_e_mpc,mpc_warm);
toc
disp('MPC:'); disp(converged_MPC);

MPC_Prob = ProbCheck(MPC_Prob,'snopt');
mpc_warm.sol = 1;
save('MPC_WARM_TMPC.mat','mpc_warm');

%% Setup CCM Controller
[geo_Prob,cntrl_info] = setup_CCM_controller(n,m,geodesic_N,...
                                          W_fnc,dW_fnc,n_W,f,B,...
                                          lambda,MP_state(1,:)',MP_ctrl(1,:)',test_state);