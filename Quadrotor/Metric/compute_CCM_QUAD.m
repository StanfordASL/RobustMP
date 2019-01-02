clear all; close all; clc;
% yalmip('clear');

warning('off','MATLAB:lang:badlyScopedReturnValue');

%% Constants

n = 9;
% lambda = 1;
g = 9.8066;
ccm_eps = 1e-1;

r_lim = pi/3;
p_lim = pi/3;
th_lim_low = 0.5*g;
th_lim_high = 2*g;
vx_lim = 1.5;
vy_lim = 1.5;
vz_lim = 1.5;

%% Uncomment for global opt

% lambda_range = linspace(0.84,0.93,3);
% lambda_range = (1/100)*round(lambda_range*100);
% euc_bounds = NaN(length(lambda_range),1);
% d_bars = NaN(length(lambda_range),1);
% cond_bound = NaN(length(lambda_range),1);
% 
% eps = 0.1;
% condn_prev = 10;
% return_metric = 0;
% 
% for ll = 1:length(lambda_range)
%     lambda = lambda_range(ll);
%     
%     fprintf('**********\n');
%     fprintf('lambda: %f\n', lambda);
%     solved = 0;
%     
%     %Determine upper bound
%     cond_l = condn_prev;
%     cond_u = 1.2*condn_prev;
%     while (~solved)
%         fprintf(' cond_u: %.2f: ', cond_u);
%         [sos_prob,~,~,~] = CCM_QUAD_Opt(n,g,r_lim,p_lim,th_lim_low,th_lim_high,vx_lim,vy_lim,vz_lim,...
%             cond_u,lambda,ccm_eps,return_metric);
%         if (sos_prob == 0)
%             solved = 1;
%             fprintf('feasible \n');
%         else
%             %shift up condition number range
%             fprintf('\n');
%             cond_l = cond_u;
%             cond_u = 1.2*cond_u;
%         end
%     end
%     if (solved)
%         euc_bounds(ll) = sqrt(cond_u)/lambda;
%         fprintf(' cond_l: %.2f, cond_u: %.2f\n', cond_l, cond_u);
%     else
%         continue;
%     end
%     
%     %Now do bisection search
%     while(cond_u - cond_l >= eps)
%         condn = (cond_l+cond_u)/2;
%         fprintf(' cond: %.2f', condn);
%         try
%             [sos_prob,w_lower,w_upper,~] = CCM_QUAD_Opt(n,g,r_lim,p_lim,th_lim_low,th_lim_high,vx_lim,vy_lim,vz_lim,...
%                 condn,lambda,ccm_eps,return_metric);
%         catch
%             sos_prob = 1;
%         end
%         if (sos_prob == 0)
%             fprintf(' feasible\n');
%             
%             euc_bounds(ll) = sqrt(condn)/lambda;
%             d_bars(ll) = sqrt(double(1/w_lower))/lambda;
%             
%             cond_u = condn;
%         else
%             fprintf(' infeasible\n');
%             cond_l = condn;
%         end
%     end
%     condn_prev = cond_u;
%     cond_bound(ll) = cond_u;
%     disp('Euc_bound:'); disp(euc_bounds(ll));
%     disp('d_bar:'); disp(d_bars(ll));
%     fprintf('**********\n');
%     
% end
% 
% save('quad_global_opt_bounds_fine.mat');

%% Pick a solution

lambda = 1.29;
condn = 150;
return_metric = 1;

save_name = 'metric_QUAD_vectorized_2.mat';
[sos_prob,w_lower,w_upper,~] = CCM_QUAD_Opt(n,g,r_lim,p_lim,th_lim_low,th_lim_high,vx_lim,vy_lim,vz_lim,...
    condn,lambda,ccm_eps,return_metric,save_name);

%% Compute aux control bound
load(save_name);
disp('Checking CCM conditions and Computing control bound...');
lambda = 0.998*lambda;

Bw = [zeros(3);
    eye(3);
    zeros(3)];
B_perp = [eye(6);
    zeros(3,6)];

ctrl_N = 15;
vx_range = linspace(-vx_lim,vx_lim,ctrl_N);
vy_range = linspace(-vy_lim,vy_lim,ctrl_N);
vz_range = linspace(-vz_lim,vz_lim,ctrl_N);
r_range = linspace(-r_lim, r_lim, ctrl_N);
p_range = linspace(-p_lim, p_lim, ctrl_N);
th_range = linspace(th_lim_low, th_lim_high, 2*ctrl_N);

b_T = @(x) [sin(x(9)); -cos(x(9))*sin(x(8)); cos(x(9))*cos(x(8))];

f_mat = @(x) [[x(4);x(5);x(6)];
    [0;0;g] - b_T(x)*x(7);
    zeros(3,1)];

%gradients
db_T_q =@(x) [0, cos(x(9));
    -cos(x(8))*cos(x(9)), sin(x(8))*sin(x(9));
    -sin(x(8))*cos(x(9)),-cos(x(8))*sin(x(9))];

%          x y z vx     vy    vz  t r p
df_perp_mat      = @(x) [zeros(3), eye(3),zeros(3,3);
    zeros(3,6), -b_T(x), -(db_T_q(x)*[1;0])*x(7), -(db_T_q(x)*[0;1])*x(7)];

%     eig_CCM = zeros(ctrl_N,ctrl_N,ctrl_N,ctrl_N, ctrl_N, ctrl_N);
%     eig_W = zeros(ctrl_N,ctrl_N, ctrl_N,ctrl_N,ctrl_N,ctrl_N,2);
%     sigma_ThBw = zeros(ctrl_N,ctrl_N,ctrl_N,ctrl_N, ctrl_N, ctrl_N);
eig_CCM = zeros(ctrl_N,ctrl_N,ctrl_N);
eig_W = zeros(ctrl_N,ctrl_N, ctrl_N,2);
sigma_ThBw = zeros(ctrl_N,ctrl_N,ctrl_N);

%     for i1 = 1:length(vx_range)
%         for i2 = 1:length(vy_range)
%             for i3 = 1:length(vz_range)
i1 = 1; i2 = 1; i3 = 1;

for i4 = 1:length(r_range)
    for i5 = 1:length(p_range)
        for i6 = 1:length(th_range)
            
            x = [zeros(3,1);
                vx_range(i1);vy_range(i2);vz_range(i3);
                th_range(i6);r_range(i4); p_range(i5)];
            
            W = W_eval(w_poly_fnc(x));
            M = W\eye(n);
            Theta = chol(M);
            Theta_Bw = Theta*Bw;
            sigma_ThBw(i4,i5,i6) = max(sqrt(eig(Theta_Bw'*Theta_Bw)));
            
            L = chol(W);
            
            f = f_mat(x);
            df_perp = df_perp_mat(x);
            
            dW_f = W_eval(dw_poly_vx_fnc(x))*f(4) + W_eval(dw_poly_vy_fnc(x))*f(5) + W_eval(dw_poly_vz_fnc(x))*f(6);
            F = -dW_f(1:6,1:6) + df_perp*W*B_perp + B_perp'*W*df_perp' + 2*lambda*W(1:6,1:6);
            
            R_CCM = -F;
            
            eig_CCM(i4,i5,i6) = min(eig(R_CCM));
            eig_W(i4,i5,i6,1) = min(eig(W));
            eig_W(i4,i5,i6,2) = max(eig(W));
        end
    end
end
%             end
%         end
%     end

d_bar = 0.1*max(sigma_ThBw(:))/lambda;
disp('d_bar'); disp(d_bar);
disp('W:'); disp(min(min(min(eig_W(:,:,:,1)))));
disp(max(max(max(eig_W(:,:,:,2)))));
disp('CCM:'); disp(min(eig_CCM(:)));

disp('euc_bounds');
disp(d_bar*sqrt(diag(W_upper_mat)));





