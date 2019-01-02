clear all; close all; clc;
% yalmip('clear');

warning('off','MATLAB:lang:badlyScopedReturnValue');

%% Constants

n = 6;
g = 9.81;
ccm_eps = 0.05;

%Dynamics constants

p_lim = pi/3;
pd_lim = pi/3;
vx_lim = 2;
vz_lim = 1.0;


%% Uncomment for global optimization of bound

% lambda_range = linspace(0.7,0.95,5);
% lambda_range = (1/100)*round(lambda_range*100);
% euc_bounds = NaN(length(lambda_range),1);
% d_bars = NaN(length(lambda_range),1);
% cond_bound = NaN(length(lambda_range),1);
% 
% eps = 1;
% condn_prev = 50;
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
%         [sos_prob,~,~] = CCM_PVTOL_Opt(n,g,p_lim,pd_lim,vx_lim,vz_lim,...
%                                 cond_u,lambda,ccm_eps,return_metric);
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
%         
%         [sos_prob, w_lower, w_upper] = CCM_PVTOL_Opt(n,g,p_lim,pd_lim,vx_lim,vz_lim,...
%                                 condn,lambda,ccm_eps,return_metric);
%         
%         if (sos_prob == 0)
%             fprintf(' feasible\n');
%             
%             euc_bounds(ll) = sqrt(double(w_upper/w_lower))/lambda;
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
% pause;

%% Pick a (lambda,condn)

lambda = 0.83; 
condn = 132.8;
return_metric = 1;

save_name = 'metric_PVTOL_vectorized.mat';
[sos_prob, w_lower, w_upper] = CCM_PVTOL_Opt(n,g,p_lim,pd_lim,vx_lim,vz_lim,...
                                condn,lambda,ccm_eps,return_metric,save_name);

%% Compute aux control bound
load(save_name);

disp('Checking CCM conditions and Computing control bound...');
lambda = 0.998*lambda;

dW_vz_fnc = @(x) zeros(n);
dW_pd_fnc = @(x) zeros(n);

m = 0.486;
J = 0.00383;
len = 0.25;

B = [zeros(1,4),1/m,len/J;
     zeros(1,4),1/m,-len/J]';
 
Bw = @(x)[zeros(1,3),cos(x(3)),-sin(x(3)),0]';
  
B_perp = [eye(4);
        zeros(2,4)];

ctrl_N = 12;
p_range = linspace(-p_lim, p_lim, ctrl_N);
vx_range = linspace(-vx_lim, vx_lim, ctrl_N);
vz_range = linspace(-vz_lim, vz_lim, ctrl_N);
pd_range = linspace(-pd_lim, pd_lim, ctrl_N);

sin_x = @(x) 0.7264*(x/(pi/4));
cos_x = @(x) 0.8516 - 0.1464*(2*(x/(pi/4))^2 -1);

df_mat = @(x) [0,0,-x(4)*sin(x(3))-x(5)*cos(x(3)),cos(x(3)),-sin(x(3)),0; 
               0,0, x(4)*cos(x(3))-x(5)*sin(x(3)),sin(x(3)), cos(x(3)),0;
               zeros(1,5),1;
               0,0,-g*cos(x(3)),0,x(6),x(5);
               0,0, g*sin(x(3)),-x(6),0,-x(4);
               zeros(1,6)];

f_mat = @(x) [x(6); %p_dot
              x(6)*x(5) - g*sin(x(3)); %vy_dot
              -x(6)*x(4) - g*cos(x(3)); %vz_dot
              0]; %pd_dot
          

delta_u = zeros(ctrl_N,ctrl_N,ctrl_N,ctrl_N);
eig_CCM = zeros(ctrl_N, ctrl_N, ctrl_N, ctrl_N);
eig_W = zeros(ctrl_N,ctrl_N,2);
sigma_ThBw = zeros(ctrl_N,ctrl_N,ctrl_N,ctrl_N);

for i = 1:length(p_range)
    for j = 1:length(vx_range)
        for k = 1:length(vz_range)
            for l = 1:length(pd_range)
                x = [randn(2,1);p_range(i);vx_range(j);vz_range(k);pd_range(l)];
                
                W = W_eval(w_poly_fnc(x));
                M = W\eye(n);
                Theta = chol(M);
                Theta_Bw = Theta*Bw(x);
                sigma_ThBw(i,j,k,l) = max(sqrt(eig(Theta_Bw'*Theta_Bw)));
                
                L = chol(W);
                
                f = f_mat(x);
                df = df_mat(x);
                F = -W_eval(dw_poly_p_fnc(x))*f(1) - W_eval(dw_poly_vy_fnc(x))*f(2)-...
                     dW_vz_fnc(x)*f(3) - dW_pd_fnc(x)*f(4) + ...
                     df*W + W*df' + 2*lambda*W;
                
                delta_u_den = eig((inv(L))'*(B*B')*inv(L));
                delta_u(i,j,k,l) = 0.5*max(eig((inv(L))'*F*inv(L)))/...
                    sqrt(min(delta_u_den(delta_u_den>0)));
                
                R_CCM = -B_perp'*F*B_perp;
                
                eig_CCM(i,j,k,l) = min(eig(R_CCM));
                eig_W(i,j,1) = min(eig(W));
                eig_W(i,j,2) = max(eig(W));
            end
        end
    end
end
alpha_w = max(sigma_ThBw(:));
d_bar = alpha_w/lambda;
disp('d_bar'); disp(d_bar);
disp('Control:'); disp(max(d_bar*delta_u(:)));
disp('W:'); disp(min(min(eig_W(:,:,1))));
disp(max(max(eig_W(:,:,2))));
disp('CCM:'); disp(min(eig_CCM(:)));

disp('euc_bounds');
disp(d_bar*sqrt(diag(W_upper)));









