function [solved, w_lower, w_upper] = CCM_Tube_Opt(x1_lim, x2_lim,condn, lambda, ccm_eps,return_metric,varargin)

norm_scale = 1e-5;

%% State-space and dynamics

x = msspoly('x',2);

f = [-1*x(1) + 2*x(2);
    -3*x(1) + 4*x(2) - 0.25*(x(2)^3)];
B = [0.5;-2];

df = [-1, 2;
     -3, 4-0.75*x(2)^2];

dtwo = msspoly('dtwo',2);

%% Initialize problem

prog = spotsosprog;
prog = prog.withIndeterminate(x);
prog = prog.withIndeterminate(dtwo);

[prog, w_lower] = prog.newPos(1);
[prog, w_upper] = prog.newPos(1);

%% Parameterize W

w_order = 4;
w_poly = monomials(x,0:w_order);
W_list = cell(length(w_poly),1);

[prog, W_list{1}] = prog.newSym(2);
W = W_list{1}*w_poly(1);

for i = 2:length(w_poly)
    [prog, W_list{i}] = prog.newSym(2);
    W = W + W_list{i}*w_poly(i);
end

[prog, W_upper] = prog.newPSD(2);

dW_f = diff(W(:),x)*f;
dW_f = reshape(dW_f,2,2);

%% Killing field conditions

dW_b = diff(W(:),x)*B; %db_dx = 0
prog = prog.withPolyEqs(dW_b);

%%

%Box constraints
box_lim = [x(1)+x1_lim;
           x1_lim-x(1);
           x(2)+x2_lim;
           x2_lim-x(2)];

%Upper and Lower definiteness
l_order = 4;
l_def_states = [x;dtwo]';

[prog, Ll] = prog.newSOSPoly(monomials(l_def_states,0:l_order),4);
[prog, Lu] = prog.newSOSPoly(monomials(l_def_states,0:l_order),4);

%CCM condition
l_ccm_states = [x;dtwo]';
lc_order = 4;

[prog, Lc] = prog.newSOSPoly(monomials(l_ccm_states,0:lc_order),4);

%W uniform bounds
% prog = prog.withPos(w_lower-0.0035);
prog = prog.withPos(w_lower-1.0);
prog = prog.withPSD(w_upper*eye(2)-W_upper);

%Condition bound
prog = prog.withPos(condn*w_lower - w_upper);
% prog = prog.withPos(condn*w_lower - trace(W_upper));

%W pos def
prog = prog.withSOS((dtwo'*W*dtwo - w_lower*(dtwo'*dtwo)) - Ll'*box_lim);
prog = prog.withSOS(dtwo'*(W_upper - W)*dtwo - Lu'*box_lim);

%CCM condition
[prog, rho] = prog.newFreePoly(monomials(x,0:4),1);

R_CCM = -(-dW_f + df*W + W*df' - rho*(B*B') + 2*lambda*W);
prog = prog.withSOS((dtwo'*R_CCM*dtwo - ccm_eps*(dtwo'*dtwo)) - Lc'*box_lim);

options = spot_sdp_default_options();
options.verbose = return_metric;

%Norm constraint
free_vars = [prog.freeVar];
len = length(free_vars);
[prog, a] = prog.newPos(len);
prog = prog.withPos(-free_vars + a);
prog = prog.withPos(free_vars + a);
try
    SOS_soln = prog.minimize(norm_scale*sum(a), @spot_mosek, options);
catch
    %failed
    solved = 1;
    w_lower = 0;
    w_upper = 0;
    return;
end
solved = ~strcmp(SOS_soln.info.solverInfo.itr.prosta, 'PRIMAL_AND_DUAL_FEASIBLE');

w_lower = double(SOS_soln.eval(w_lower));
w_upper = double(SOS_soln.eval(w_upper));

if (return_metric)
    if (solved==0)
        disp('feasible, getting results...');
        
        W_sol = zeros(2,2,length(w_poly));
        NNZ_list = zeros(length(w_poly),1);
        for i = 1:length(w_poly)
            W_sol(:,:,i) = clean(double(SOS_soln.eval(W_list{i})),1e-4);
            if sum(sum(abs(W_sol(:,:,i)))) > 0
                NNZ_list(i) = 1;
            end
        end
        w_poly = w_poly(find(NNZ_list));
        W_sol = W_sol(:,:,find(NNZ_list));
        fprintf('%d unique monomials\n',length(w_poly));
        
        dw_poly_x1 = diff(w_poly,x(1));
        dw_poly_x2 = diff(w_poly,x(2));
        
        W_upper = clean(double(SOS_soln.eval(W_upper)),1e-4);
        
%         pause;
        
        %% Create monomial functions
        w_poly_fnc = mss2fnc(w_poly,x,randn(length(x),2));
        dw_poly_x1_fnc = mss2fnc(dw_poly_x1,x,randn(length(x),2));
        dw_poly_x2_fnc = mss2fnc(dw_poly_x2,x,randn(length(x),2));
        
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
        save(varargin{1},'W_eval','w_poly_fnc','dw_poly_x1_fnc','dw_poly_x2_fnc','W_upper');
        
    end
end