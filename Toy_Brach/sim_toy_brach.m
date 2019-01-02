clear; close all; clc;

%% Toy Brach Problem

n = 2; %state dimen
m = 1; %control dim
N = 10;

%Number of collocation points
K = N;

%CGL nodes
[s_t,w] = clencurt(K); %t_t: [-1, 1] : <-> : [0, Tp]
s = fliplr(s_t); %t: [1, -1]


%% Get Differentiation matrix

D = ChebyshevDiffMatrix(N,s); %arranged for forward time
D = kron(D,eye(n));

%% Variables

%State node values: [x_t0,...,x_tN]
%Control node values: [u_t0,...,u_tN]

% n_vars = (N+1)*(n+m);

%% Setup problem

Aeq = [eye(n), zeros(n,n*N+N+1+1);
       zeros(1,n*N),1,0,zeros(1,N+1+1)];
   
b_L = [zeros(2,1);0.5];
b_U = b_L;  

nl_constr = @(xt) brach_con(xt,D,N,n);
% nl_conJ = @(xt) brach_cond(xt,D,N,n);
%n_vars = 3*(N+1)+1;

x0 = [kron([zeros(N,1);0.5],[1;0])+...
      kron((linspace(0,1,N+1)+0.1)',[0;1]);
      zeros(N+1,1);0.4];
x_L = [zeros(2*(N+1),1);-pi*ones(N+1,1);0];
x_U = [5*ones(2*(N+1),1);pi*ones(N+1,1);2];

Name = 'Brach';
Prob = conAssign(@(xt) xt(end),[],[],[],...
            x_L,x_U,Name, x0,...
            [], 0, Aeq, b_L,b_U,...
            nl_constr, [],[],[],...
            zeros(2*(N+1),1),zeros(2*(N+1),1),...
            [],[],[],[]);
        
%% Solve

Prob = ProbCheck(Prob,'snopt');

Result = snoptTL(Prob);

converged = Result.Inform; %GOOD: {1,2,3}

%% Interpret solution

%Final time
Tp = Result.x_k(end);

%evaluate solution at grid tau
delta = Tp;
dt = 0.01;
tau = 0:dt:Tp;
s_e = (2*tau - Tp)/Tp; %[-1, s_delta]

%CGL time points
tau_CGL = (Tp*s_t+Tp)/2;

%Lagrange polynomial evaluation at the interpolation points
L_e = compute_Lagrange(length(s_e)-1,N,s_e,s_t);

%Compute trajectories
NMPC_state = zeros(size(L_e,2),n);
x_nom = zeros(N+1,n);
for i = 1:n
    c = Result.x_k(i:n:n*(N+1)-(n-i))';
    NMPC_state(:,i) = (c*L_e)';
    x_nom(:,i) = c';
end


NMPC_ctrl = zeros(size(L_e,2),m);
u_nom = zeros(N+1,m);
for j = 1:m
    c = Result.x_k(n*(N+1)+j:m:end-1-(m-j))';
    NMPC_ctrl(:,j) = (c*L_e)';
    u_nom(:,j) = c';
end

%% Plot

figure()
plot(tau',NMPC_state,'linewidth',2);
hold on
plot(tau_CGL,x_nom,'rs','markersize',10);
grid on
xlabel('time [s]');
legend('x','y');

figure()
plot(tau',NMPC_ctrl,'linewidth',2);
hold on
plot(tau_CGL,u_nom,'rs','markersize',10);
grid on
xlabel('time [s]');
ylabel('\theta');