clear all; close all; clc;
yalmip('clear');

%% Constants

x1_lim = 5;
x2_lim = 5;
ccm_eps = 0.01;

%% Uncomment for global optimization

% lambda_range = linspace(0.8,2,15); %line search range
% euc_bounds = NaN(length(lambda_range),1);
% d_bars = NaN(length(lambda_range),1);
% cond_bound = NaN(length(lambda_range),1);
% 
% eps = 1e-3;
% condn_prev = 1.63;
% return_metric = 0;

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
%         [sos_prob,~,~] = CCM_Tube_Opt(x1_lim, x2_lim,cond_u, lambda, ccm_eps,return_metric);
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
%         fprintf(' cond_l: %.4f, cond_u: %.4f\n', cond_l, cond_u);
%     else
%         continue;
%     end
%     
%     %Now do bisection search
%     while(cond_u - cond_l >= eps)
%         condn = (cond_l+cond_u)/2;
%         fprintf(' cond: %.4f', condn);
%         
%         [sos_prob, w_lower, w_upper] = CCM_Tube_Opt(x1_lim, x2_lim,condn, lambda, ccm_eps,return_metric);
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

% figure()
% plot(lambda_range, euc_bounds,'ro','markerfacecolor','g','markersize',20);
% grid on
% xlabel('\lambda');
% ylabel('$\|x^{*}-x\|/\bar{w}$','Interpreter','Latex');
% % title('Robustness optimization');
% set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)

% pause;

%% Pick a solution

lambda = 1.742857142857143; 
condn = 1.6344;

% lambda = 1.8;
% condn = 3.1;
return_metric = 1;

save_name = 'metric_tube.mat';
[sos_prob, w_lower, w_upper] = CCM_Tube_Opt(x1_lim, x2_lim,condn, lambda, ccm_eps,return_metric,save_name);

%% Compute control bounds using optimal metric (W = const matrix)
load(save_name);

f_mat = @(x) [-1*x(1) + 2*x(2);
    -3*x(1) + 4*x(2) - 0.25*(x(2)^3)];

df_mat = @(x) [-1, 2;
     -3, 4-0.75*x(2)^2];
B = [0.5;-2];
B_perp = [2;0.5];
Bw = [0;1];

x1_range = linspace(-x1_lim,x1_lim,30);
x2_range = linspace(-x2_lim,x2_lim,30);
delta_u = zeros(length(x2_range),length(x1_range));
eig_CCM = delta_u;
sigma_ThBw = zeros(length(x2_range),length(x1_range));

for i = 1:length(x2_range)
    for j = 1:length(x1_range)
        x = [x1_range(j); x2_range(i)];
        
        W = W_eval(w_poly_fnc(x));
        M = W\eye(2);
        sigma_ThBw(i,j) = max(sqrt(eig(Bw'*M*Bw)));
        
        L = chol(W);
        f = f_mat(x);
        df = df_mat(x);
       
        F_sol = -W_eval(dw_poly_x1_fnc(x))*f(1) - W_eval(dw_poly_x2_fnc(x))*f(2) + ...
            df_mat(x)*W + W*df_mat(x)' + 2*lambda*W;
    
        delta_u(i,j) = 0.5*max(eig((inv(L))'*F_sol*inv(L)))/...
                          sqrt(max(eig((inv(L))'*(B*B')*inv(L))));
                      
        eig_CCM(i,j) = min(eig(B_perp'*(-F_sol)*B_perp));
    end
end

w = 0.1;
d_bar = max(sigma_ThBw(:))*w/lambda;

figure()
surf(x1_range,x2_range,d_bar*delta_u);
grid on
xlabel('x1'); ylabel('x2'); zlabel('$\bar{\delta}_u$','interpreter','latex');
set(findall(gcf,'type','text'),'FontSize',28);set(gca,'FontSize',28)


disp('d_bar:'); disp(d_bar);
disp('control:'); disp(max(delta_u(:))*d_bar);
disp('CCM:'); disp(min(eig_CCM(:)));

if length(w_poly_fnc([1;1])) == 1
    M = W_eval(w_poly_fnc([1;1]))\eye(2); %W= const. otherwise - replace with W_upper
else
    M = W_upper\eye(2);
end

%Compare RCI sets
P_rci = diag([39.0251, 486.0402]);

figure()
Ellipse_plot((W_upper\eye(2))*(1/d_bar^2),[0;0],20,'c'); hold on
Ellipse_plot(M*(1/d_bar^2),[0;0],20,'k',0.6); 
Ellipse_plot(P_rci,[0;0],20,'r');

