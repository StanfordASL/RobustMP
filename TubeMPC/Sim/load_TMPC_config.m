%% Constants
n = 2;
m = 1;

%% Setup Metric

load '../Metric/metric_tube.mat';

w_poly = @(x) w_poly_fnc(x);
W_fnc = struct('W_eval',W_eval,'w_poly_fnc',w_poly);
dW_fnc = @(x) {}; %expect W to be constant

sigma_ThBw = 0.355;
lambda =  1.7429;
ctrl_bound = 1.908;
n_W = []; %states that W is a function of

%% Dynamics & Cost

f = @(x) [-1*x(1) + 2*x(2); 
          -3*x(1) + 4*x(2) - 0.25*(x(2)^3)];
df = @(x) [-1, 2;
           -3, 4-0.75*(x(2)^2)];
B = @(x) [0.5;-2];

dB = @(x) zeros(n,n,m);

B_w = [0;1];

Q = diag([0.5;0.5]); R = 1;

%% Bounds

w_max = 0.1;

%expecting W to be a constant - else replace with W_upper\eye(n)
M_ccm = W_fnc.W_eval(1)\eye(n); 
d_bar = (w_max*sigma_ThBw/lambda);
ctrl_bound = ctrl_bound*w_max;
euc_bound = d_bar*sqrt(diag(W_fnc.W_eval(1)));

%terminal set
alpha = 10;
P = [7.9997, -12.2019;
    -12.2019, 27.0777];

%% Simulation constraints

state_constr_low = -[5;5]+euc_bound;
state_constr = [state_constr_low, -state_constr_low];
ctrl_constr = [-2+ctrl_bound,2-ctrl_bound];

x_eq = [0;0];
u_eq = 0;
test_state = [3.4;-2.4];

