%% Constants

Jx = 0.03; Jy =  0.04; Jz =  0.1;
global mq Jq;

Jq = diag([Jx;Jy;Jz]);
mq = 0.9574;

g = 9.81;

n = 12;
m = 4;

%% Dynamics

%full-state dynamics
%pos,vel, roll, pitch, yaw, om
f = @(x) [x(4:6);
          [0;0;g];
          R_eul(x(7:9))*x(10:12);
          -Jq\cross(x(10:12),Jq*x(10:12))];
      
B =  @(x) [zeros(3,4);
          -rot_matrix(x(7),x(8),x(9))*[0;0;1],zeros(3);
          zeros(3,4);
          zeros(3,1), Jq\eye(3)];  

B_w = [zeros(3);
       eye(3);
       zeros(6,3)];
  
%dynamics for CCM controller
%pos,vel,thrust,roll,pitch
f_ctrl = @ f_xc;      
B_ctrl = @(xc)[zeros(6,4);
               eye(4)];     

%% Setup metric & bounds

%disturbance bound
w_max = 0.1; 

%metric type
pullback = 0;

%Metric Controller state-space: xc: x,y,z, vx,vy,vz, th, r,p,y
%Metric Controller control-space: uc: th_dot, rd, pd, yd
nc = 10; mc = 4;
 
load ('../Metric/metric_QUAD_vectorized_2.mat');

W_fnc = struct('W_eval',W_eval,'w_poly_fnc',w_poly_fnc);
dW_fnc = @(x) {dw_poly_th_fnc(x), dw_poly_r_fnc(x), dw_poly_p_fnc(x)};
n_W = [7,8,9]; %states W is a function of

sigma_ThBw = 0.5564;
lambda = 1.29;

%Bounds
W_upper = W_upper_mat;
d_bar = (w_max*sigma_ThBw/lambda);

In = eye(9);
%Maximal position tube
M_ccm_pos = (1/d_bar^2)*((In(1:3,:)*W_upper*In(1:3,:)')\eye(3));
[U_pos,S_pos,V_pos] = svd(M_ccm_pos);
S_new = (sqrt(S_pos\eye(3)) + 0.2*eye(3))^2\eye(3);
M_ccm_pos_infl = U_pos*S_new*V_pos';

%choose a bound for yaw
yaw_bound = 10*(pi/180);
%choose metric component for yaw
M_ccm_yaw = (d_bar/yaw_bound)^2;
W_ccm_yaw = (1/M_ccm_yaw);
   
%% Setup lower-level controller

%Angular rate P controller gain
global kp_om;
kp_om = 6*lambda;

%% Setup Sim environment
close all;

%Set world space dimensions
World_dim = [5,5,5];

%set initial state
pos_init = [0;0;1]; %in NWU frame
test_state = [pos_init(1);-pos_init(2);-pos_init(3); %NED frame
              zeros(9,1)];
thrust_init = g;          
xc_init = [test_state(1:6);
           thrust_init;
           test_state(7:9)];
          
%Define obstacle environment (if not already done) in NWU frame
if (new_obs)
    define_obstacles;
    save('create_env/Obstacles.mat','tree_obstacles','tower_obstacles',...
          'obstacles','obstacles_infl','obstacles_coll','n_obstacles');
else
    load('create_env/Obstacles.mat');
    fig = figure();
    plot_all_obstacles();
end  

%define goal
goal_box = [World_dim-[1,1,2];
            World_dim];
goal_V =  corner2Hull(goal_box,1);
goal = Polyhedron('V',goal_V);
figure(fig)
plot3dObstacles(goal_box,'b',0.2);

