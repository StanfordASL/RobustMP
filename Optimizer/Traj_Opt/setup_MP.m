function [MP_Prob,L_e_full] = ...
    setup_MP(n,m,...
               f,B,df,dB,x_con,u_con,...
               N,Tp,dt,...
               P,alpha,RPI_bound,...
               Q,x_eq,R,u_eq,obs)

%Create global motion planning problem (usually needs a decent seed)

%%%% Inputs %%%%
%n,m: state and control dimension
%f,B,df,dB: f (n x 1),B (n x m),dfdx (n x n),dBdx (n x n x m) function handles
%x_con,u_con: lower and upper state and control bounds 
%N: #collocation points
%Tp: time-horizon (potentially can make a variable)
%dt: time-resolution of final solution
%x_eq,u_eq: final state&control values (think equilibrium)
%P,alpha: terminal constraint: \|x*(Tp) - x_eq\|_P^2 \leq \alpha
%RPI_bound: initial state bound based on geodesic energy:  E(x*(0),x(0)) <= RPI_bound
%Q,R: state and control cost matrices
%obs: obstacle information struct (if relevant)

%%%% Outputs %%%%
%MP_prob: tomlab problem struct
%L_e_full: Lagrange interpolating polynomial for final solution

%% Constants

%State bounds
x_L = x_con(:,1);
x_U = x_con(:,2);

%Control bounds
u_L = u_con(:,1);
u_U = u_con(:,2);

%Number of collocation points-1
K = N;

%CGL nodes
[s_t,w] = clencurt(K); %t_t: [-1, 1] : <-> : [0, Tp]
s = fliplr(s_t); %t: [1, -1]

%% Final solution interpolation matrix

tau_full = 0:dt:Tp; 
s_e_full = (2*tau_full - Tp)/Tp; %[-1, 1]

%Lagrange polynomial evaluation at the interpolation points (if want to
%define additional constraints along trajectory)
% L_e = compute_Lagrange(length(s_e)-1,N,s_e,s_t);

%Evaluation at full grid (finer than collocation grid)
L_e_full = compute_Lagrange(length(s_e_full)-1,N,s_e_full,s_t);

%% Get Differentiation matrix

D = ChebyshevDiffMatrix(N,s); %arranged for forward time
D = kron(D,eye(n));

%% Variables

%State node values: [x_t0,...,x_tN]
%Control node values: [u_t0,...,u_tN]

% n_vars = (N+1)*(n+m);

%% Define problem

x_eq_all = kron(ones(N+1,1),x_eq);
u_eq_all = kron(ones(N+1,1),u_eq);

xu_eq = [x_eq_all;u_eq_all];

Q_bar = kron(diag(w),Q); R_bar = kron(diag(w),R);
Q_tilde = Q_bar + kron(diag([zeros(N,1);(2/Tp)]),P);

F = blkdiag(Q_tilde,R_bar);
% F_pattern = sparse(F~=0);
    
xu_L = [kron(ones(N+1,1),x_L);
        kron(ones(N+1,1),u_L)];
xu_U = [kron(ones(N+1,1),x_U);
        kron(ones(N+1,1),u_U)]; 
       
MP_cost = @(xu) (Tp/2)*(xu-xu_eq)'*F*(xu-xu_eq);
MP_grad = @(xu) Tp*F*(xu-xu_eq);
MP_hess = @(xu) Tp*F;

%constraints:
%dynamics, initial RPI constraint, terminal set constraint, obstacles
c_L = [zeros(n*(N+1)+2,1); ones(obs.n_obs*(N+1),1)];
c_U = [zeros(n*(N+1),1);RPI_bound;alpha;Inf*ones(obs.n_obs*(N+1),1)];

con = @(xu,Prob) MP_con(xu,Prob,n,m,N,P,D,f,B,Tp,x_eq,obs);
conJ = @(xu,Prob) MP_conJ(xu,Prob,n,m,N,P,D,df,dB,Tp,x_eq,obs);

%initial guess placeholder
xu0 = zeros((n+m)*(N+1),1);

Name = 'MP';
MP_Prob = conAssign(MP_cost,MP_grad,MP_hess,[],...
            xu_L,xu_U,Name, xu0,...
            [], 0, [],[],[],...
            con,conJ,[],[],...
            c_L,c_U,...
            [],[],[],[]);
        
        
% NMPC_Prob.SOL.optPar(10) = 1e-3;
        
MP_Prob.user.x_act = zeros(n,1);

end