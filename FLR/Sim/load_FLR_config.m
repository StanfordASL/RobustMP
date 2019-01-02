%% Constants

n = 4;
m = 1;

%% Dynamics 

mass = 1;
l = 1;
I = (1/3)*mass*(2*l)^2;
sigma = 100;
J = 1;
b = 1;
g = 9.81;

f  = @(x) [x(2);
          (mass*g*l/I)*sin(x(1)) - (sigma/I)*(x(1)-x(3));
           x(4);
           (sigma/J)*(x(1)-x(3)) - (b/J)*x(4)];
       
B = @(x) [zeros(3,1);1/J];

df = @(x) [0, 1, 0, 0;
           (mass*g*l/I)*cos(x(1))-(sigma/I), 0, (sigma/I), 0;
           zeros(1,3), 1;
           (sigma/J), 0, -(sigma/J), -(b/J)];
       
dB = @(x) zeros(n,n,m);       

B_w = [0, (1/I), 0, 0;
       0, 0, 0, (1/J)]';

%% Obstacle info

obs = struct('n_obs',0); %there are no obstacles

%% Setup Metric & constants

load '../Metric/metric_FLR_vectorized.mat';

w_poly = @(x) w_poly_fnc(wrapToPi(x(1,:)));
W_fnc = struct('W_eval',W_eval,'w_poly_fnc',w_poly);
dW_fnc = @(x) {dw_poly_x1_fnc(wrapToPi(x(1,:)))};

sigma_ThBw = 0.0475;
lambda =  2.5;
ctrl_bound = 5.91; %normalized
n_W = 1; %states that W is a function of
   
%% Setup Bounds

w_max = 1;

d_bar = (w_max*sigma_ThBw/lambda);
ctrl_bound = ctrl_bound*w_max;
euc_bound = d_bar*sqrt(diag(W_upper));

In = eye(n);
M_ccm_ang = (1/d_bar^2)*((In([1;3],:)*W_upper*In([1;3],:)')\eye(2));
M_ccm_omg = (1/d_bar^2)*((In([2;4],:)*W_upper*In([2;4],:)')\eye(2));

state_constr_low = -[pi;5;pi;5]+euc_bound;
x_con = [state_constr_low, -state_constr_low];
u_con = [-35+ctrl_bound, 35-ctrl_bound];

%% Setup sim environment
           
q_eq = 0*(pi/180); %straight up
th_eq = q_eq - (mass*g*l/sigma)*sin(q_eq);
x_eq = [q_eq; 0; th_eq; 0];
u_eq = -sigma*(q_eq-th_eq); 

link_ang = -pi+15*(pi/180);
mot_ang = link_ang + 5*(pi/180);
test_state = [link_ang;
              -30*(pi/180);
              mot_ang;
              0];