function [LMP_Prob,L_e_full,s_t] = ...
    setup_LMP(MP_x,MP_u,MP_t,n,m,...
               f,B,df,dB,x_con,u_con,...
               N,T,delta,dt,...
               P,alpha,RPI_bound,...
               Q,x_eq,R,u_eq,mu,obs)
%Create local receding-horizon motion planning problem (LMP)

%%%% Inputs %%%%
%MP_x, MP_u, MP_t: global motion plan state, control, time trajectory
%n,m: state and control dimension
%f,B,df,dB: f (n x 1),B (n x m),dfdx (n x n),dBdx (n x n x m) function handles
%x_con,u_con: lower and upper state and control bounds 
%N: #collocation points
%T: MPC time-horizon 
%delta: re-solve time for MPC
%dt: time-resolution of final solution
%x_eq,u_eq: final state&control values (think equilibrium)
%P,alpha: terminal constraint: \|x*(Tp) - MP_x(T*)\|_P^2 \leq \alpha
%RPI_bound: initial state bound based on geodesic energy:  E(x*(0),x(0)) <= RPI_bound
%Q,R: state and control cost matrices
%mu: weight for re-join time
%obs: obstacle information struct (if relevant)

%%%% Outputs %%%%
%NMPC_prob: tomlab problem struct
%L_e_full: Lagrange interpolating polynomial for final solution over [0,T]
%s_t: normalized [-1,1] time nodes for collocation points
%% Constants

%State bounds
x_L = x_con(:,1);
x_U = x_con(:,2);

%Control bounds
u_L = u_con(:,1);
u_U = u_con(:,2);

%Number of collocation points
K = N;

%CGL nodes
[s_t,w] = clencurt(K); %s_t: [-1, 1] : <-> : [0, Tp]
s = fliplr(s_t); %t: [1, -1]

%% Final solution interpolation matrix

tau_full = 0:dt:T; 
s_e_full = (2*tau_full - T)/T; %[-1, 1]

%Lagrange polynomial evaluation at the interpolation points
L_e_full = compute_Lagrange(length(s_e_full)-1,N,s_e_full,s_t);

%% Get Differentiation matrix

D = ChebyshevDiffMatrix(N,s); %arranged for forward time
D = sparse(kron(D,eye(n)));

%% Variables

%State node values: [x_t0,...,x_tN]
%Control node values: [u_t0,...,u_tN]
%Rejoin MP time: T_star

% n_vars = (N+1)*(n+m) + 1;

%% Define problem

% u_eq = zeros(m,1);

x_eq_all = kron(ones(N+1,1),x_eq);
u_eq_all = kron(ones(N+1,1),u_eq);

xu_eq = [x_eq_all;u_eq_all];

Q_bar = kron(diag(w),Q); R_bar = kron(diag(w),R);
F = sparse(blkdiag(Q_bar,R_bar));

xu_L = [kron(ones(N+1,1),x_L);
        kron(ones(N+1,1),u_L);
        delta];
xu_U = [kron(ones(N+1,1),x_U);
        kron(ones(N+1,1),u_U);
        MP_t(end)]; 
       
MPC_cost = @(xu) (T/2)*(xu(1:(N+1)*(n+m))-xu_eq)'*F*(xu(1:(N+1)*(n+m))-xu_eq)-mu*xu(end);
MPC_grad = @(xu) [T*F*(xu(1:(N+1)*(n+m))-xu_eq);-mu];
MPC_hess = @(xu) blkdiag(T*F,0);

%constraints:
%dynamics,initial,terminal,obstacles
c_L = [zeros(n*(N+1)+2,1); ones(obs.n_obs*(N+1),1)];
c_U = [zeros(n*(N+1),1);RPI_bound;alpha;Inf*ones(obs.n_obs*(N+1),1)];

MPC_con = @(xu,Prob) LMP_con(xu,Prob,n,m,N,P,D,f,B,T,MP_t,MP_x);
MPC_conJ = @(xu,Prob) LMP_conJ(xu,Prob,n,m,N,P,D,f,df,B,dB,T,MP_t,MP_x,MP_u);

%initial guess placeholder
xu0 = [zeros((n+m)*(N+1),1);T];

global LMP_CON_J; %initialize for efficiency
LMP_CON_J = zeros(n*(N+1)+2+obs.n_obs*(N+1),(n+m)*(N+1)+1);

Name = 'LMP';
LMP_Prob = conAssign(MPC_cost,MPC_grad,MPC_hess,[],...
            xu_L,xu_U,Name, xu0,...
            [], 0, [],[],[],...
            MPC_con,MPC_conJ,[],[],...
            c_L,c_U,...
            [],[],[],[]);
        
        
LMP_Prob.SOL.optPar(10) = 1e-4;
LMP_Prob.SOL.optPar(12) = 1e-4;
 
LMP_Prob.user.x_act = zeros(n,1);

end