function [geo_Prob,T_e,T_dot_e] = ...
    setup_geodesic_calc(n,N,W,dW,n_W)

%%%%% Inputs %%%%%
%n: state-space dimension
%N: Chebyshev polynomial order
%W, dW: W and dW matrix functions
%n_W: states that W is a function of

%%%%% Outputs %%%%%
%geo_Prob: Tomlab problem structure
%K_e: #evaluation points for final solution
%T_e: evaluation matrix for geodesic curve
%T_dot_e: evaluation matrix for geodesic curve velocity


%% Define constants

K = 2*N;
K_e = 3; %beginning,middle,endpoint

%Optimization variables: chebyshev coefficients for geodesic
%{c_0^1...c_N^1},...,{c_0^1...c_N^n}

%% Obtain Chebyschev Pseudospectral Numerics

%CGL points and quadrature weights
[t,w] = clencurt(K);
[t_e,~] = clencurt(K_e-1);

%compute Cheby basis at all points
[T, T_dot] = ...
    compute_cheby(K,N,t);

[T_e, T_dot_e] = ...
    compute_cheby(K_e-1,N,t_e);

%% Geodesic endpoint constraints

[phi_start, ~] = compute_cheby(0,N,-1);
[phi_end, ~] = compute_cheby(0,N,1);
A_start = kron(eye(n),phi_start');
A_end = kron(eye(n),phi_end');

Aeq = sparse([A_start;
              A_end]);

%% For evaluating x_k, x_dot_k

Phi = zeros(n,n*(N+1),K+1);
Phi_dot = zeros(n,n*(N+1),K+1);
for k = 1:K+1
    Phi(:,:,k) = kron(eye(n),T(:,k)');
    Phi_dot(:,:,k) = 2*kron(eye(n),T_dot(:,k)');
end

%use to evaluate cost derivative
Ti = zeros(n*(N+1),K+1,length(n_W));
In = eye(n);
for j = 1:length(n_W)
    i = n_W(j);
    Ti(:,:,j) = kron(In(:,i),T);
end

%% Setup cost and gradient functions

global  GEO_X; 
GEO_X = zeros(n,K+1);

global  GEO_MXDOT; 
GEO_MXDOT = zeros(n,K+1);

% Cost function
geo_cost_fnc =  @(vars) Geodesic_cost(vars,w,n,...
    K,W,Phi,Phi_dot);

%Gradient function
geo_grad_fnc = @(vars) Geodesic_grad(vars,w,...
    K,N,n,Ti,W,dW,Phi,Phi_dot,n_W);

%% Define problem

Name = 'Geodesic Problem';
geo_Prob = conAssign(geo_cost_fnc,geo_grad_fnc,[],[],...
                  [],[],Name,zeros(n*(N+1),1),[],0,...
                  Aeq,zeros(2*n,1),zeros(2*n,1),[],[],[],[],[],[]);

end