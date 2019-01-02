clear all; close all; clc;
warning('off','MATLAB:nargchk:deprecated');
         
%% Load constants, dynamics, CCM, bounds, sim initial condns

load_FLR_config;

%% Load & initialize solvers

%load and initialize MP and CCM problem structs
load_solvers;

%% Visualize

visualize_FLR;

%% Set up non-linear sim

ode_options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);

dt_sim = 0.002;
t_end = Tp;
solve_t = (0:dt_sim:t_end)';
T_steps = length(solve_t)-1;

w_dist = zeros(T_steps,2);

X = zeros(T_steps+1,n);

U_fb = zeros(T_steps,m);
U = zeros((t_end/dt)+1,m);

%computation time
ctrl_solve_time = NaN(T_steps,1);

%solve success
opt_solved = NaN(T_steps,1);

%geodesic energy
geo_energy = zeros(T_steps,1);

X(1,:) = test_state';
x = test_state;

%% Simulate
disp('Ready to Simulate');
keyboard;

for i = 1:T_steps
    
    if (mod(i,500)==0)
        fprintf('%d / %d\n',i,T_steps);
    end
    
    x_nom = MP_state(1+(i-1)*(dt_sim/dt),:);
    u_nom = MP_ctrl(1+(i-1)*(dt_sim/dt):1+i*(dt_sim/dt),:); 
    
    %Feedback control
    tic
    [J_opt,opt_solved(i),geo_Prob,cntrl_info,u_fb] = compute_CCM_controller(geo_Prob,cntrl_info,...
        x_nom',u_nom(1,:)',x);
    ctrl_solve_time(i,1) = toc;
    
    geo_energy(i) = J_opt;
    
    U_fb(i,:) = u_fb';
    
    U(1+(i-1)*(dt_sim/dt):1+i*(dt_sim/dt),:) = u_nom+repmat(U_fb(i,:),(dt_sim/dt)+1,1);
    
    %Disturbance model
    w_dist(i,:) = -(w_max/sqrt(2))*[1,1];
    
    %sim
    [d_t,d_state] = ode113(@(t,d_state)ode_sim(t,d_state,[solve_t(i):dt:solve_t(i+1)]',u_nom,U_fb(i,:),...
        f,B,B_w,w_dist(i,:)'),[solve_t(i),solve_t(i+1)],x,ode_options);
    
    %update
    x = d_state(end,:)';
    X(i+1,:) = x';
end

%% Plots
close all;
plot_FLR;

%% Movie sim

keyboard;
plot_FLR_movie(solve_t,X,round(0.05/dt_sim));
