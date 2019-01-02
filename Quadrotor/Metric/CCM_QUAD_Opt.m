function [solved,w_lower,w_upper,W_upper_mat] = ...
    CCM_QUAD_Opt(n,g,r_lim,p_lim,th_lim_low,th_lim_high,vx_lim,vy_lim,vz_lim,...
            condn,lambda,ccm_eps,return_metric,varargin)
%%


% W_scale = (1e-2)*diag([0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.0001;0.01;0.01]);
W_scale = (1e-5)*diag([12;12;24;0.1;0.1;0.2;5;3.5;3.5]);
% W_scale = zeros(9);
norm_scale = 1e-7;

% sin_x = @(x) 0.5059*(x/(pi/6));
% cos_x = @(x) 0.9326 - 0.06699*(2*(x/(pi/6))^2 -1);

sin_x = @(x) 0.9101*(x/(pi/3)) - 0.04466*(4*(x/(pi/3))^3 - 3*(x/(pi/3)));
cos_x = @(x) 0.7441 -0.2499*(2*(x/(pi/3))^2 -1);

%states
n = 9;
x = msspoly('x',9);

%pos_def indeterminates
dnin = msspoly('dnin',9);
dsix = msspoly('dsix',6);

sin_r = sin_x(x(8));
cos_r = cos_x(x(8));
sin_p = sin_x(x(9));
cos_p = cos_x(x(9));

%dynamics f
b_T = [sin_p; -cos_p*sin_r; cos_p*cos_r];

f = [[x(4);x(5);x(6)];
     [0;0;g] - b_T*x(7);
     zeros(3,1)]; 

%gradients       
db_T_q = [0, cos_p;
         -cos_r*cos_p, sin_r*sin_p;
         -sin_r*cos_p,-cos_r*sin_p];

%          x y z vx     vy    vz t r p     
df_perp      = [zeros(3), eye(3),zeros(3,3);
                zeros(3,6),-b_T, -db_T_q(:,1)*x(7), -db_T_q(:,2)*x(7)];
            
B_perp = [eye(6);
          zeros(3,6)];
      
%% Initialize problem

prog = spotsosprog;
prog = prog.withIndeterminate(x);
prog = prog.withIndeterminate(dnin);
prog = prog.withIndeterminate(dsix);

[prog, w_lower] = prog.newPos(1);
[prog, w_upper] = prog.newPos(1);
[prog, W_upper] = prog.newSym(n);

%% Parametrize W (2)

w_states = x(7:9);

w_order = 4;
[w_poly, w_poly_mat] = monomials(w_states,0:w_order);

W_perp_list = cell(length(w_poly),1);
W_list = cell(length(w_poly),1);
W_pc_list = cell(length(w_poly),1);
W_c_list = cell(length(w_poly),1);

[prog, W_perp_list{1}] = prog.newSym(6);
[prog, W_pc_list{1}] = prog.newFree(6,3);
[prog, W_c_list{1}] = prog.newSym(3);

W_list{1} = [W_perp_list{1}, W_pc_list{1};
            W_pc_list{1}', W_c_list{1}];
W = W_list{1}*w_poly(1);

for i = 2:length(w_poly)
    if size(w_poly_mat,2)>3 %more than (r,p,th)
        if (sum(w_poly_mat(i,end-2:end))>0)
            W_perp_list{i} = zeros(6);
        else
            [prog, W_perp_list{i}] = prog.newSym(6);
        end
    else
        W_perp_list{i} = zeros(6);
    end
    [prog, W_c_list{i}] = prog.newSym(3);
    [prog, W_pc_list{i}] = prog.newFree(6,3);
    
    W_list{i} = [W_perp_list{i}, W_pc_list{i};
                 W_pc_list{i}', W_c_list{i}];
    W = W + W_list{i}*w_poly(i);
end

W_perp = W(1:6,1:6);
dW_perp_f = diff(W_perp(:),x)*f;
dW_perp_f = reshape(dW_perp_f,6,6);

%% Definiteness conditions

%Lagrange multipliers
box_lim = [r_lim^2-x(8)^2;
           p_lim^2-x(9)^2;
           x(7) - th_lim_low;
           th_lim_high - x(7)];
%            vz_lim^2-x(6)^2];

l_order = w_order;
l_def_states = w_states;
n_def_L = length(box_lim);

[w_def_mon, w_def_mat] = monomials([l_def_states;dnin],0:l_order);
w_def_keep = find(sum(w_def_mat(:,length(l_def_states)+1:end),2)==2); %only keep quadratics in dnin
w_def_mon = w_def_mon(w_def_keep);

[prog, Ll] = prog.newSOSPoly(w_def_mon,n_def_L);
[prog, Lu] = prog.newSOSPoly(w_def_mon,n_def_L);

l_ccm_states = w_states;
lc_order = l_order;

[ccm_def_mon_rp, ccm_def_mat_rp] = monomials([l_ccm_states;dsix],0:lc_order+2);
ccm_def_keep_rp = find(sum(ccm_def_mat_rp(:,length(l_ccm_states)+1:end),2)==2); %only keep quadratics in dsix
ccm_def_mon_rp = ccm_def_mon_rp(ccm_def_keep_rp);

[ccm_def_mon_th, ccm_def_mat_th] = monomials([l_ccm_states;dsix],0:lc_order+2);
ccm_def_keep_th = find(sum(ccm_def_mat_th(:,length(l_ccm_states)+1:end),2)==2); %only keep quadratics in dsix
ccm_def_mon_th = ccm_def_mon_th(ccm_def_keep_th);

% [prog, Lc_v]  = prog.newSDSOSPoly(monomials(l_ccm_states,0:2),3);
[prog, Lc_rp] = prog.newSOSPoly(ccm_def_mon_rp,2);
[prog, Lc_th] = prog.newSOSPoly(ccm_def_mon_th,length(box_lim)-2);
Lc = [Lc_rp; Lc_th];

%W uniform bounds
prog = prog.withPos(w_lower-1);
prog = prog.withPSD(w_upper*eye(n)-W_upper);

%Condition bound
prog = prog.withPos(condn*w_lower - w_upper);

%W pos def
prog = prog.withSOS((dnin'*W*dnin - w_lower*(dnin'*dnin)) - (Ll'*box_lim(1:n_def_L)));
prog = prog.withSOS(dnin'*(W_upper - W)*dnin - (Lu'*box_lim(1:n_def_L)));

%CCM condition
R_CCM = -(-dW_perp_f + df_perp*W*B_perp + B_perp'*W*df_perp' + 2*lambda*W_perp);
prog = prog.withSOS((dsix'*R_CCM*dsix - ccm_eps*(dsix'*dsix)) - (Lc'*box_lim));

options = spot_sdp_default_options();
options.verbose = return_metric;

%Norm constraint
free_vars = [prog.coneVar(2:end); prog.freeVar];
len = length(free_vars);
[prog, a] = prog.newPos(len);
prog = prog.withPos(-free_vars + a);
prog = prog.withPos(free_vars + a);

try
    SOS_soln = prog.minimize(norm_scale*sum(a) + trace(W_scale*W_upper), @spot_mosek, options);
catch
    %failed
    solved = 1;
    w_lower = 0;
    w_upper = 0;
    W_upper_mat = zeros(n);
    return;
end

try
    solved = ~(strcmp(SOS_soln.info.solverInfo.itr.prosta, 'PRIMAL_AND_DUAL_FEASIBLE'));% && ...
%            strcmp(SOS_soln.info.solverInfo.itr.solsta, 'OPTIMAL'));
catch
    solved = 1;
    w_lower = 0;
    w_upper = 0;
    W_upper_mat = zeros(n);
    return;
end

%% Parse

if (solved == 0)
    w_lower = double(SOS_soln.eval(w_lower));
    w_upper = double(SOS_soln.eval(w_upper));
    W_upper_mat = clean(double(SOS_soln.eval(W_upper)),1e-4);
else
    w_lower = 0;
    w_upper = 0; 
    W_upper_mat = zeros(n);
    return;
end

if (return_metric)
    if (solved==0)
        disp('feasible, getting results...');

        W_sol = zeros(n,n,length(w_poly));
        NNZ_list = zeros(length(w_poly),1);
        for i = 1:length(w_poly)
            W_sol(:,:,i) = clean(double(SOS_soln.eval(W_list{i})),1e-7);
            if sum(sum(abs(W_sol(:,:,i)))) > 0
                NNZ_list(i) = 1;
            end
        end
        w_poly = w_poly(find(NNZ_list));
        W_sol = W_sol(:,:,find(NNZ_list));
        
        fprintf('%d non-zero monomials\n',length(w_poly));
        
        %save in coefficient form (for copying over to C++)
        p = w_poly_mat(find(NNZ_list),:);
        save('Quad_Metric_Coeffs.mat','W_sol','p');
        
        dw_poly_th = diff(w_poly,x(7));
        dw_poly_r = diff(w_poly,x(8));
        dw_poly_p = diff(w_poly,x(9));
        dw_poly_vx = diff(w_poly,x(4));
        dw_poly_vy = diff(w_poly,x(5));
        dw_poly_vz = diff(w_poly,x(6));
        
%         pause;
        
        %% Create monomial functions
        w_poly_fnc = mss2fnc(w_poly,x,randn(length(x),2));
        dw_poly_r_fnc = mss2fnc(dw_poly_r,x,randn(length(x),2));
        dw_poly_p_fnc = mss2fnc(dw_poly_p,x,randn(length(x),2));
        dw_poly_th_fnc = mss2fnc(dw_poly_th,x,randn(length(x),2));
        dw_poly_vx_fnc = mss2fnc(dw_poly_vx,x,randn(length(x),2));
        dw_poly_vy_fnc = mss2fnc(dw_poly_vy,x,randn(length(x),2));
        dw_poly_vz_fnc = mss2fnc(dw_poly_vz,x,randn(length(x),2));
        
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
        save(varargin{1},'W_eval','w_poly_fnc','dw_poly_r_fnc','dw_poly_p_fnc','dw_poly_th_fnc','dw_poly_vx_fnc','dw_poly_vy_fnc','dw_poly_vz_fnc','W_upper_mat');
    end
end
end