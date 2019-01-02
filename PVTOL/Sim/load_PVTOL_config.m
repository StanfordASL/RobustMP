%% Dimensions

n = 6;
m = 2;

%% Dynamics

mass = 0.486;
J = 0.00383;
g = 9.81;
len = 0.25;

f  = @f_pvtol;
       
B = @(x)[zeros(1,4),1/mass, len/J;
         zeros(1,4),1/mass,-len/J]';

df = @df_pvtol;

dB = @(x) zeros(n,n,m);

B_w = [zeros(1,3),1,0,0;
       zeros(1,3),0,1,0]';

%% Setup Metric & constants

load '../Metric/metric_PVTOL_vectorized.mat';

W_fnc = struct('W_eval',W_eval,'w_poly_fnc',w_poly_fnc);
dW_fnc = @(x) {dw_poly_p_fnc(x), dw_poly_vy_fnc(x)};

alpha_w = 0.3127;
lambda =  0.8283;
ctrl_bound = 6.00; %normalized
n_W = [3,4]; %states that W is a fnc of

%% Setup Bounds

%disturbance bound
w_max = 0.1;

%geodesic distance bound
d_bar = (w_max*alpha_w/lambda);

%control and euclidean bounds
ctrl_bound = ctrl_bound*w_max;
euc_bound = d_bar*sqrt(diag(W_upper));

%state and control constraint adjustment
state_constr_low = -[5.5;5.5;pi/4;2;1;pi/3]+euc_bound;
x_con = [state_constr_low, -state_constr_low];
u_con = [0.1*mass*g+ctrl_bound, 2*mass*g-ctrl_bound;
         0.1*mass*g+ctrl_bound, 2*mass*g-ctrl_bound];

%% Setup sim environment

%target state     
x_eq = [4.5;4.5;0;0;0;0];

%equilibium control
u_eq = [0.5*mass*g; 0.5*mass*g]; 

%initial state
test_state = [-4.4;
              -5;
               0;
               0.5;
               0;
               0];

%all obstacles
obs_loc = [[-4;-3.5],...
           [3;-4],...
           [0.7;-3],...
          [-1;-0.5],...
           [2.5;-0.5],...
           [-4;1],...
           [1.0;1.7],...
           [2.5;3.8],...
           [-2;4]];

obs_rad = [0.5,1,0.9,0.8,0.9,1,0.9,0.5,0.6];       

obs = struct('n_obs',length(obs_rad),'pos',obs_loc,'r',obs_rad);

In = eye(n);
%Unscaled position tube
M_ccm_pos_unscaled = ((In(1:2,:)*W_upper*In(1:2,:)')\eye(2)); 
%Maximal position tube
M_ccm_pos = (1/d_bar^2)*((In(1:2,:)*W_upper*In(1:2,:)')\eye(2)); 
[U_pos,S_pos,V_pos] = svd(M_ccm_pos);

%Rescale maximal position tube by obstacle + robot radius
M_obs = zeros(2,2,obs.n_obs);
for i = 1:obs.n_obs
    S_new = (sqrt(S_pos\eye(2)) + (obs_rad(i)+len)*eye(2))^2\eye(2);
    M_obs(:,:,i) = U_pos*S_new*V_pos';
end
obs.M_obs = M_obs;
