function [MPC_state,MPC_ctrl,converged,warm,Prob] = ...
    compute_NMPC(Prob,x_act,x_con,u_con,x_eq,u_eq,...
    n,m,N,L_e,warm)
%Solve receding horizon MPC problem

%%%% Inputs %%%%%
%Prob: tomlab problem struct
%x_act: current state
%mpc_state: guess for starting MP state (usually just equal to x_act)
%x_con: state constraints
%u_con: control constraints
%x_eq,u_eq: terminal equilibrium state and control
%n,m: state/control dimensions
%N: number of collocation points for solution
%L_e: Lagrange interpolating polynomial for solution
%warm: warm-start struct

%%%% Outputs %%%%%
%MPC_state,MPC_ctrl: state/ctrl plan
%converged: flag
%warm: warm-start solution for future use
%Prob: updated problem struct for future use

%% Solution guess
if (~warm.sol) %don't have warm-start guess
    %this should only happen for the first MPC iteration
    %use straight line guess
    X0 = zeros(N+1,n);
    U0 = repmat(u_eq',N+1,1);
    for i = 1:n
        X0(:,i) = linspace(x_act(i),x_eq(i),N+1);
    end
    x0 = reshape(X0',(N+1)*n,1);
    u0 = reshape(U0',(N+1)*m,1);
    
    Prob = modify_x_0(Prob,[x0;u0]);
    
else
    %recall warm soln
    Prob = WarmDefSOL('snopt',Prob,warm.result);
end

%% Update constraint information

%update initial state
Prob.user.x_act = x_act; 

%% Pre-check

if ~Prob.CHECK
    Prob = ProbCheck(Prob,'snopt');
end

%% Solve
Result = snoptTL(Prob);

converged = Result.Inform; %GOOD: {1,2,3}

%% Compute trajectories
MPC_state = zeros(size(L_e,2),n);
x_nom = zeros(N+1,n);
for i = 1:n
    c = Result.x_k(i:n:n*(N+1)-(n-i))';
    MPC_state(:,i) = (c*L_e)';
    x_nom(:,i) = c';
end


MPC_ctrl = zeros(size(L_e,2),m);
u_nom = zeros(N+1,m);
for j = 1:m
    c = Result.x_k(n*(N+1)+j:m:end-(m-j))';
    MPC_ctrl(:,j) = (c*L_e)';
    u_nom(:,j) = c';
end

if (converged==1 || converged==2 || converged==3) %only save converged solutions as warm starts
    warm.result = Result;
    warm.state = x_nom;
    warm.ctrl = u_nom;
    warm.sol = 1;
else
    warm.sol = 0;
end


end