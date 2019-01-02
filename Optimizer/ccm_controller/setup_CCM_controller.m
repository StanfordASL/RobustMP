function [geo_Prob,cntrl_info] = setup_CCM_controller(n,m,geodesic_N,...
                                                      W_fun,dW_fun,n_W,f_fun,B_fun,...
                                                      lambda,x_0,u_0,x)

% setup CCM controller 

%%%%% Inputs %%%%%
%n,m: state and control dimension
%geodesic_N: order of Chebyshev approximation for geodesic
%W_fun, dW_fun: function handles for W(x) and dW_dx
%n_W: variables W is a function of
%f_fun,B_fun function handles for f(x), B(x)
%lambda: contraction rate
%x_0, u_0: test initial nominal state and control
%x: test initial state

%%%%% Outputs %%%%%
% geo_Prob: geodesic problem structure (Tomlab);
% cntrl_inf: controller struct

%% Setup geodesic solver

[geo_Prob,T_e,T_dot_e] = setup_geodesic_calc(n,geodesic_N,W_fun,dW_fun,n_W);

%choose solver
geo_solver = 'npsol';   

%initialize solution structure (used later for warm-start)
geo_warm = struct('sol',0,'result',[],'E',0);

%% define controller struct 
cntrl_info.n = n;
cntrl_info.m = m;
cntrl_info.T_e = T_e;
cntrl_info.T_dot_e = T_dot_e;
cntrl_info.N = geodesic_N;
cntrl_info.solver = geo_solver;
cntrl_info.warm = geo_warm;
cntrl_info.lambda = lambda;
cntrl_info.W = W_fun;
cntrl_info.dW = dW_fun;
cntrl_info.f = f_fun;
cntrl_info.B = B_fun;

%% First time run (to generate solution structure)

% start_p = n x 1 start point of geodesic (nominal state)
% end_p = n x 1 end point of geodesic (actual state)

tic
[J_opt,converged_geo,geo_Prob,cntrl_info,u_aux]  = compute_CCM_controller(geo_Prob,cntrl_info,x_0,u_0,x);
toc;

fprintf('Geodesic:%d (good: 0, 1, 6), Geodesic dist: %.3f, u_fb: \n',converged_geo,sqrt(J_opt));
disp(u_aux);
geo_Prob.CHECK = 1;

end