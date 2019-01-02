function [MP_state,MP_ctrl,converged,warm] = ...
    compute_MP(Prob,x_act,x_con,u_con,x_eq,u_eq,...
                 n,m,N,L_e,warm)

%Solve global motion planning problem

%%%% Inputs %%%%%
%Prob: tomlab problem struct
%x_act: current state
%x_con: state constraints
%u_con: control constraints
%x_eq,u_eq: terminal equilibrium state and control
%n,m: state/control dimensions
%N: number of collocation points for solution
%L_e: Lagrange interpolating polynomial for solution
%warm: warm-start struct

%%%% Outputs %%%%%
%MP_state,MP_ctrl: state/ctrl plan
%converged: flag
%warm: warm-start solution for future use

%% Guess or recall

if (~warm.sol) %Solution guess
    mp_state = x_act;
    for i = 1:n
        if mp_state(i) > x_con(i,2)
            mp_state(i) = x_con(i,2)-0.01*(x_con(i,2)-state_const(i,1));
        elseif mp_state(i) < x_con(i,1)
            mp_state(i) = x_con(i,1)+0.01*(x_con(i,2)-state_const(i,1));
        end
    end
    In = eye(n);
    x0 = zeros((N+1)*n,1);
    for i = 1:n
        x0 = x0 + kron(linspace(mp_state(i),x_eq(i),N+1), In(i,:))';
    end
    
    Im = eye(m);
    u0 = zeros((N+1)*m,1);
    for j = 1:m
        u0 = u0 + kron(u_eq(j)*ones(N+1,1), Im(:,j));
    end
    
    Prob = modify_x_0(Prob,[x0;u0]);
else %recall warm solution
    Prob = WarmDefSOL('snopt',Prob,warm.result);
end

%% Update constraint information
Prob.user.x_act = x_act;
Prob = ProbCheck(Prob,'snopt');

%% Solve
Result = snoptTL(Prob);

converged = Result.Inform; %GOOD: {1,2,3}

%% Compute trajectories
MP_state = zeros(size(L_e,2),n);
x_nom = zeros(N+1,n);
for i = 1:n
    c = Result.x_k(i:n:n*(N+1)-(n-i))';
    MP_state(:,i) = (c*L_e)';
    x_nom(:,i) = c';
end


MP_ctrl = zeros(size(L_e,2),m);
u_nom = zeros(N+1,m);
for j = 1:m
    c = Result.x_k(n*(N+1)+j:m:end-(m-j))';
    MP_ctrl(:,j) = (c*L_e)';
    u_nom(:,j) = c';
end

warm.result = Result;
warm.state = x_nom;
warm.ctrl = u_nom;

end