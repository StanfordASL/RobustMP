function [MPC_state,MPC_ctrl,T_rejoin,converged,warm,Prob,obs_active] = ...
    compute_LMP(Prob,x_act,x_con,u_con,MP_t,MP_x,MP_u,...
    n,m,N,L_e,warm,dt,t_s,delta,E_s,lambda,obs)
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
%dt: time-resolution of solution
%t_s: solve time
%delta: MPC re-solve time
%E_s: updated upper bound on on E(x*(0),x(0))
%lambda: contraction rate
%obs: full obstacle set

%%%% Outputs %%%%%
%MPC_state,MPC_ctrl: state/ctrl plan
%T_rejoin: optimal rejoin time T*
%converged: flag
%warm: warm-start solution for future use
%Prob: updated problem struct for future use
%obs_active: active set of obstacles considered

%% Solution guess
if (~warm.sol) %don't have warm-start guess
    %this should only happen for the first MPC iteration
    %use nominal MP
    tau = min(t_s + (1/2)*(warm.Tp*warm.s_t+warm.Tp),MP_t(end)) ;
    state_nom = interp1(MP_t,MP_x,tau);
    ctrl_nom = interp1(MP_t,MP_u,tau);
    
    x0 = zeros((N+1)*n,1);
    u0 = zeros((N+1)*m,1);
    for k = 1:N+1
        x_nom = state_nom(k,:)';
        for i = 1:n
            if x_nom(i) >= x_con(i,2)
                x_nom(i) = x_con(i,2)-0.001*(x_con(i,2)-x_con(i,1));
            elseif x_nom(i) < x_con(i,1)
                x_nom(i) = x_con(i,1)+0.001*(x_con(i,2)-x_con(i,1));
            end
        end
        
        x0(1+(k-1)*n:k*n) = x_nom;
        
        u_nom = ctrl_nom(k,:)';
        for j = 1:m
            if (u_nom(j) <= u_con(j,1))
                u_nom(j) = u_con(j,1)+0.001*(u_con(j,2)-u_con(j,1));
            elseif (u_nom(j) >= u_con(j,2))
                u_nom(j) = u_con(j,2)-0.001*(u_con(j,2)-u_con(j,1));
            end
        end
        u0(1+(k-1)*m:k*m) = u_nom;
    end

%     x_guess = zeros(n,N+1);
%     for i = 1:n
%         x_guess(i,:) = linspace(x_act(i),state_nom(end,i),N+1);
%     end
%     u_guess = repmat(mean(ctrl_nom)',1,N+1);
%     
%     x0 = reshape(x_guess,(N+1)*n,1);
%     u0 = reshape(u_guess,(N+1)*m,1);
    
    Prob = modify_x_0(Prob,[x0;u0;warm.Tp]);
    
else
    %recall warm soln
    Prob = WarmDefSOL('snopt',Prob,warm.result);
end

%% Update constraint information

%update initial state
Prob.user.x_act = x_act; 

%update initial bound
Prob = modify_c_U(Prob,E_s,n*(N+1)+1);

%lower bound on nominal rejoin time = solve_time + delta
Prob = modify_x_L(Prob,min(t_s+delta,MP_t(end)),(n+m)*(N+1)+1);

%% Update obstacle set (prune and dynamic re-size)

if obs.n_obs > 0
    obs_mpc = obs;
    
    %First let's prune the obstacles we need to look at
    t_ind_l = (t_s/dt)+1;
    t_ind_u = min((t_s+1.5*warm.Tp)/dt,size(MP_x,1)); %conservative upper-bound on T_rejoin
    
    x_min = min(MP_x(t_ind_l:t_ind_u,1:2));
    x_max = max(MP_x(t_ind_l:t_ind_u,1:2));
    box_min = x_min'-max(obs_mpc.r)*[1;1];
    box_max = x_max'+max(obs_mpc.r)*[1;1];
    
    obs_active = zeros(obs_mpc.n_obs,1);
    for i = 1:obs_mpc.n_obs
        if ( sum(obs_mpc.pos(:,i)>=box_min)==2 && ...
                sum(obs_mpc.pos(:,i)<=box_max)==2 )
            %obstacle within range we care about
            obs_active(i) = 1;
            
            if obs_mpc.time_var
                %do time re-scale of obstacles
                tau = (1/2)*(warm.Tp*warm.s_t+warm.Tp);
                E_time_bound = (sqrt(E_s)*exp(-lambda*tau) + warm.d_bar*(1-exp(-lambda*tau))).^2;
                for k = 1:N+1
                    S_new = (sqrt(E_time_bound(k)*(obs_mpc.S\eye(2))) + (obs_mpc.r(i))*eye(2))^2\eye(2);
                    obs_mpc.M_obs(:,:,i,k) = obs_mpc.U*S_new*obs_mpc.V';
                end
            end
%         else
%             %set its location far away so the constraint becomes irrelevant
%             obs_mpc.pos(:,i) = [100;100];
        end
    end
    
    Prob.user.obs = obs_mpc;
else
    obs_active = [];
end

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
    c = Result.x_k(n*(N+1)+j:m:end-1-(m-j))';
    MPC_ctrl(:,j) = (c*L_e)';
    u_nom(:,j) = c';
end

T_rejoin = Result.x_k(end);

warm.result = Result;
warm.state = x_nom;
warm.ctrl = u_nom;
warm.sol = 1;

end