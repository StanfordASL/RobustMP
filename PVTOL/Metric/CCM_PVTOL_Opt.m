function [solved,w_lower,w_upper] = ...
    CCM_PVTOL_Opt(n,g,p_lim,pd_lim,vx_lim,vz_lim,...
    condn,lambda,ccm_eps,return_metric,varargin)
%%


% W_scale = diag([0.02;0.04;0.06;0.005;0.01;0.008]);
W_scale = diag([0.001;0.02;0.0367;0.005;0.001;0.002]);
norm_scale = 1e-4;

% sin_x = @(x) 0.05059*(x/(pi/6));
% cos_x = @(x) 0.9326 - 0.06699*(2*(x/(pi/6))^2 -1);

% sin_x = @(x)  0.7264*(x/(pi/4)) - 0.01942*(4*(x/(pi/4))^3 - 3*(x/(pi/4)));
% cos_x = @(x) 0.8516 - 0.1464*(2*(x/(pi/4))^2 -1);

sin_x = @(x) 0.9101*(x/(pi/3)) - 0.04466*(4*(x/(pi/3))^3 - 3*(x/(pi/3)));
cos_x = @(x) 0.7441 -0.2499*(2*(x/(pi/3))^2 -1);

%states
x = msspoly('x',6);

%pos_def indeterminates
dsix = msspoly('dsix',6);
dfor = msspoly('dfor',4);

sin_p = sin_x(x(3));
cos_p = cos_x(x(3));

%dynamics f
f = [x(4)*cos_p - x(5)*sin_p;
    x(4)*sin_p + x(5)*cos_p;
    x(6);
    x(6)*x(5)-g*sin_p;
    -x(6)*x(4)-g*cos_p;
    0];

%          X Z  p                  x(4)    x(5)   x(6)
df_perp = [0,0,-x(4)*sin_p-x(5)*cos_p,cos_p,-sin_p,0;
           0,0, x(4)*cos_p-x(5)*sin_p,sin_p,cos_p,0;
           0,0,0,0,0,1;
           0,0,-g*cos_p,0,x(6),x(5)];

B_perp = [eye(4);
    zeros(2,4)];

%% Initialize problem

prog = spotsosprog;
prog = prog.withIndeterminate(x);
prog = prog.withIndeterminate(dsix);
prog = prog.withIndeterminate(dfor);

[prog, w_lower] = prog.newPos(1);
[prog, w_upper] = prog.newPos(1);

%% Parametrize W (2)

w_states = x(3:4);
w_poly = monomials(w_states,0:4);

W_list = cell(length(w_poly),1);
W_perp_list = cell(length(w_poly),1);
W_pc_list = cell(length(w_poly),1);
W_c_list = cell(length(w_poly),1);

[prog, W_perp_list{1}] = prog.newSym(4);
[prog, W_pc_list{1}] = prog.newFree(4,2);
[prog, W_c_list{1}] = prog.newSym(2);
W_list{1} = [W_perp_list{1},W_pc_list{1};
             W_pc_list{1}', W_c_list{1}];

W = W_list{1}*w_poly(1);

for i = 2:length(w_poly)
    [prog, W_perp_list{i}] = prog.newSym(4);
    [prog, W_pc_list{i}] = prog.newFree(4,2);
    [prog, W_c_list{i}] = prog.newSym(2);
    W_list{i} = [W_perp_list{i},W_pc_list{i};
        W_pc_list{i}', W_c_list{i}];
    
    W = W + W_list{i}*w_poly(i);
end

W_perp = W(1:4,1:4);
dW_f_perp = diff(W_perp(:),x)*f;
dW_f_perp = reshape(dW_f_perp,4,4);

[prog, W_upper] = prog.newSym(n);

%% Definiteness conditions

%Lagrange multipliers
box_lim = [ p_lim^2-x(3)^2;
            vx_lim^2-x(4)^2;
            vz_lim^2-x(5)^2;
            pd_lim^2-x(6)^2];

l_order = 4;

[pos_def_mon, pos_def_mat] = monomials([w_states;dsix],0:l_order);
pos_def_keep = find(sum(pos_def_mat(:,length(w_states)+1:end),2)==2); %only keep quadratics in dsix
pos_def_mon = pos_def_mon(pos_def_keep);
[prog, Ll] = prog.newSOSPoly(pos_def_mon,2);
[prog, Lu] = prog.newSOSPoly(pos_def_mon,2);

lc_order = 6;
l_ccm_states = [x(3);x(4);x(6);x(5)];

[ccm_def_mon, ccm_def_mat] = monomials([l_ccm_states;dfor],0:lc_order);
ccm_def_keep = find(sum(ccm_def_mat(:,length(l_ccm_states)+1:end),2)==2); %only keep quadratics in dfor
ccm_def_mon = ccm_def_mon(ccm_def_keep);
[prog, Lc] = prog.newSOSPoly(ccm_def_mon,4);

%W uniform bounds
prog = prog.withPos(w_lower-1);
prog = prog.withPSD(w_upper*eye(n)-W_upper);

%Condition bound
prog = prog.withPos(condn*w_lower - w_upper);

%W pos def
prog = prog.withSOS((dsix'*W*dsix - w_lower*(dsix'*dsix)) - (Ll'*box_lim(1:2)));
prog = prog.withSOS(dsix'*(W_upper - W)*dsix - (Lu'*box_lim(1:2)));

%CCM condition
R_CCM = -(-dW_f_perp + df_perp*W*B_perp + B_perp'*W*df_perp' + 2*lambda*W_perp);
prog = prog.withSOS((dfor'*R_CCM*dfor - ccm_eps*(dfor'*dfor)) - (Lc'*box_lim));

options = spot_sdp_default_options();
options.verbose = return_metric;

%Norm constraint
free_vars = [prog.coneVar(2:end); prog.freeVar];
len = length(free_vars);
[prog, a] = prog.newPos(len);
prog = prog.withPos(-free_vars + a);
prog = prog.withPos(free_vars + a);

try
    SOS_soln = prog.minimize(trace(W_scale*W_upper) + norm_scale*sum(a), @spot_mosek, options);
catch
    %failed
    solved = 1;
    w_lower = 0;
    w_upper = 0;
    return;
end

solved = ~SOS_soln.status.strcmp('STATUS_PRIMAL_AND_DUAL_FEASIBLE');

%% Parse

w_lower = double(SOS_soln.eval(w_lower));
w_upper = double(SOS_soln.eval(w_upper));

W_upper_mat = zeros(n);

if (return_metric)
    if (solved==0)
        disp('feasible, getting results...');

        W_sol = zeros(n,n,length(w_poly));
        NNZ_list = zeros(length(w_poly),1);
        for i = 1:length(w_poly)
            W_sol(:,:,i) = clean(double(SOS_soln.eval(W_list{i})),1e-3);
            if sum(sum(abs(W_sol(:,:,i)))) > 0
                NNZ_list(i) = 1;
            end
        end
        w_poly = w_poly(find(NNZ_list));
        W_sol = W_sol(:,:,find(NNZ_list));
        
        fprintf('%d non-zero monomials\n',length(w_poly));
        
        
        dw_poly_p = diff(w_poly,x(3));
        dw_poly_vy = diff(w_poly,x(4));
        
        W_upper = clean(double(SOS_soln.eval(W_upper)),1e-3);
        
%         pause;
        
        %% Create monomial functions
        w_poly_fnc = mss2fnc(w_poly,x,randn(length(x),2));
        dw_poly_p_fnc = mss2fnc(dw_poly_p,x,randn(length(x),2));
        dw_poly_vy_fnc = mss2fnc(dw_poly_vy,x,randn(length(x),2));
        
        %% Put together
        W_exec = 'W_eval = @(ml)';
        
        for i = 1:length(w_poly)
            if i<length(w_poly)
                W_exec = strcat(W_exec,sprintf('W_sol(:,:,%d)*ml(%d) +',i,i));
            else
                W_exec = strcat(W_exec,sprintf('W_sol(:,:,%d)*ml(%d);',i,i));
            end
        end

        %% Execute
        eval(W_exec);
        save(varargin{1},'W_eval','w_poly_fnc','dw_poly_p_fnc','dw_poly_vy_fnc','W_upper');
    end
end
end