function [J_opt,converged,Prob,controller,u_aux] = compute_CCM_controller(geo_Prob,info,x_nom,u_nom,x)

%%%%% Inputs %%%%%%
%geo_Prob: geodesic problem struct
%info: controller struct
%(x_nom,u_nom): nominal state, control
% x: current actual state

%%%%% Outputs %%%%%
%J_opt: geodesic energy
%converged: 0,1,6: good
%Prob: return controller struct with updated soln
%u_aux: feedback control

%% Compute geodesic

[X_e,X_dot_e,J_opt,converged,geo_result,Prob] = ...
    compute_geodesic(geo_Prob,info.n,info.N,x_nom,x,info.T_e,info.T_dot_e,info.warm,info.solver);

%% Update struct for next warm start

% fprintf('ok: %d, E: %.4f\n',converged, J_opt);

controller = info;
if (converged == 0 || converged == 1 || converged == 6) 
    controller.warm.sol = 1;
    controller.warm.result = geo_result;
    controller.warm.E = J_opt;
else
    disp('warn:geodesic solver failed!');
    controller.warm.sol = 0;
    u_aux = zeros(info.m,1);
    return;
end

%% Compute control

%only if successful geodesic solve 
u_aux = compute_opt_aux(X_e,X_dot_e,J_opt,info.W,info.f,info.B,u_nom,info.lambda);

end