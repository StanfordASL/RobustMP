clear; close all; clc;
warning('off','MATLAB:nargchk:deprecated');

addpath('Quad_Traj_Opt');
addpath('Quad_aux_files');
addpath('create_env');
addpath('quad_plotting');

%% Load constants, dynamics, CCM, bounds, setup sim environment

%set to 1 to generate new obs
new_obs = 0;
load_quad_config;

%% Load solver/trajectory generator

%set to 1 to  generate new path and trajectory (set 1 if new_obs = 1)
gen_new_traj = new_obs;
load_solvers;

%% Setup non-linear sim

ode_options = odeset('RelTol', 1e-6, 'AbsTol', 1e-9);

%zero-order-hold for feedback controller
dt_sim = 1/250;
t_end = t(end);
solve_t = (0:dt_sim:t_end)';
T_steps = length(solve_t)-1;

%Store disturbance
w_dist = zeros(T_steps,3);

%Store actual state
X = zeros((solve_t(end)/dt)+1,n);
X_steps = T_steps*(dt_sim/dt)+1;
Xc = zeros(T_steps+1,10);

%store errors
err = zeros(T_steps,4);

%Store control history
U_fb = zeros(T_steps,m);
U = zeros(size(X,1),m);

%Computation times
ctrl_solve_time = NaN(T_steps);

%Solve success
opt_solved = NaN(T_steps);

%Geodesic distances
geo_energy = NaN(T_steps);

%Pre-gen disturbance sequence
wind_polar = 45;
wind_az = 135;
direction = [sind(wind_polar)*cosd(wind_az),...
             sind(wind_polar)*sind(wind_az),...
             cosd(wind_polar)];
w_dist = repmat(w_max*direction,T_steps,1);


%Initialize
X(1,:) = test_state';
x = test_state;
Xc(1,:) = xc_init';
xc = xc_init;

%% Simulate

disp('Ready to Simulate');
keyboard;

for i = 1:T_steps
    if (mod(i,100)==0)
        fprintf('%d/%d \n',i,T_steps);
    end
    
    xc_nom = MP_state(1+(i-1)*(dt_sim/dt),:);
    uc_nom = MP_ctrl(1+(i-1)*(dt_sim/dt):1+i*(dt_sim/dt),:);
    
    err(i,1:3) = xc_nom(1:3)-xc(1:3)';
    err(i,4) = xc_nom(10)-xc(10);
    
    %Feedback Control
    
    tic
    [J_opt,opt_solved(i),geo_Prob,cntrl_info,u_fb] = compute_Quad_CCM_controller(geo_Prob,cntrl_info,...
        xc_nom',uc_nom(1,:)',xc);
    ctrl_solve_time(i) = toc;
    
    geo_energy(i) = J_opt;
    
    U_fb(i,:) = u_fb';
    U(1+(i-1)*(dt_sim/dt):1+i*(dt_sim/dt),:) = uc_nom+repmat(U_fb(i,:),(dt_sim/dt)+1,1);
    
    %update thrust
    xc(7) = xc(7)+(uc_nom(1,1)+u_fb(1))*dt_sim;
    
    %propagate
    [d_t,d_state] = ode113(@(t,d_state)quad_ode(t,d_state,[solve_t(i):dt:solve_t(i+1)]',uc_nom(:,2:4),...
        u_fb(2:4)',xc(7),f,B,B_w,w_dist(i,:)'),[solve_t(i):dt:solve_t(i+1)],x,ode_options);
    
    %update and record
    x = d_state(end,:)';
    X(1+(i-1)*(dt_sim/dt):1+i*(dt_sim/dt),:) = d_state;
    xc = [x(1:6);
        xc(7);
        x(7:9)];
    Xc(i+1,:) = xc';
end


%% Plot

close all
plot_quad;


