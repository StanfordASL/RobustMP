function [solved, w_lower, w_upper] = CCM_FLR_Opt(n,m,g,l,I,J,b,sigma,x1_lim,x2_lim,x3_lim,...
    condn, lambda, ccm_eps,return_metric,varargin)

W_scale = diag([0.01; 0.0005; 0.01; 0.0005]);
% W_scale = diag(zeros(4,1));

sin_x = @(x) 0.5692*(x/pi) - 0.6669*(4*(x/pi)^3 - 3*(x/pi)) + ...
    0.1043*(16*(x/pi)^5 - 20*(x/pi)^3 + 5*(x/pi));

cos_x = @(x)   -0.3042 - 0.9709*(2*(x/pi)^2 -1) + 0.3028*(8*(x/pi)^4 - 8*(x/pi)^2 + 1) + ...
    -0.029*(32*(x/pi)^6 - 48*(x/pi)^4 + 18*(x/pi)^2 - 1);

%% State-space and dynamics

%states
x = msspoly('x',2);

%pos_def indeterminates
dtre = msspoly('dtre',3);
dfor = msspoly('dfor',4);

sin_x1 = sin_x(x(1));
cos_x1 = cos_x(x(1));

%dynamics f (only need f(1) since W = W(x1))
f = x(2);
%     (m*g*l/I)*sin_x1 - (sigma/I)*(x(1)-x(3));
%     x(4);
%     (sigma/J)*(x(1)-x(3)) - (b/J)*x(4)];

df_perp = [0, 1, 0, 0;
    (m*g*l/I)*cos_x1-(sigma/I), 0, (sigma/I), 0;
    zeros(1,3), 1];

B_perp = [eye(3);
    zeros(1,3)];

%% Initialize problem

prog = spotsosprog;
prog = prog.withIndeterminate(x);
prog = prog.withIndeterminate(dtre);
prog = prog.withIndeterminate(dfor);

[prog, w_lower] = prog.newPos(1);
[prog, w_upper] = prog.newPos(1);

%% Parametrize W 

w_states = x(1);

w_order = 2;
w_poly = monomials(w_states,0:w_order);
W_list = cell(length(w_poly),1);

[prog, W_list{1}] = prog.newSym(n);
W = W_list{1}*w_poly(1);

for i = 2:length(w_poly)
    [prog, W_list{i}] = prog.newSym(n);
    W = W + W_list{i}*w_poly(i);
end

[prog, W_upper] = prog.newPSD(n);

dW_f = diff(W(:),x(1))*f;
dW_f = reshape(dW_f,n,n);

%%

%Lagrange multipliers
box_lim = [x(1)+x1_lim; x1_lim-x(1);
           x(2)+x2_lim; x2_lim-x(2)];
%            x3+x3_lim; x3_lim-x3];

l_order = 4;
l_def_states = [w_states'; dfor]';

[prog, Ll] = prog.newSOSPoly(monomials(l_def_states,0:l_order),2);
[prog, Lu] = prog.newSOSPoly(monomials(l_def_states,0:l_order),2);

l_ccm_states = [x(1);x(2);dtre]';
lc_order = 6;

[prog, Lc] = prog.newSOSPoly(monomials(l_ccm_states,0:lc_order),4);

%W uniform bounds
prog = prog.withPos(w_lower-1);
prog = prog.withPSD(w_upper*eye(n)-W_upper);

%Condition bound
prog = prog.withPos(condn*w_lower - w_upper);

%W pos def
prog = prog.withSOS((dfor'*W*dfor - w_lower*(dfor'*dfor)) - Ll'*box_lim(1:2));
prog = prog.withSOS(dfor'*(W_upper - W)*dfor - Lu'*box_lim(1:2));

%CCM condition
R_CCM = -(-dW_f(1:3,1:3) + df_perp*W*B_perp + B_perp'*W*df_perp' + 2*lambda*W(1:3,1:3));
prog = prog.withSOS((dtre'*R_CCM*dtre - ccm_eps*(dtre'*dtre)) - Lc'*box_lim);

options = spot_sdp_default_options();
% options.solveroptions.MSK_IPAR_BI_CLEAN_OPTIMIZER = 'MSK_OPTIMIZER_INTPNT';
% options.solveroptions.MSK_IPAR_INTPNT_BASIS = 'MSK_BI_NEVER';
options.verbose = return_metric;

%Norm constraint
free_vars = [prog.coneVar; prog.freeVar];
len = length(free_vars);
[prog, a] = prog.newPos(len);
prog = prog.withPos(-free_vars + a);
prog = prog.withPos(free_vars + a);

SOS_soln = prog.minimize(trace(W_scale*W_upper) + (1e-3)*sum(a), @spot_mosek, options);
try
    solved = ~SOS_soln.status.strcmp('STATUS_PRIMAL_AND_DUAL_FEASIBLE');
catch
    %failed
    solved = 1;
    w_lower = 0;
    w_upper = 0;
    return;
end

w_lower = double(SOS_soln.eval(w_lower));
w_upper = double(SOS_soln.eval(w_upper));

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
        
        dw_poly_x1 = diff(w_poly,x(1));
        
        W_upper = clean(double(SOS_soln.eval(W_upper)),1e-3);

%         pause;
        
        %% Create monomial functions
        w_poly_fnc = mss2fnc(w_poly,x,randn(length(x),2));
        dw_poly_x1_fnc = mss2fnc(dw_poly_x1,x,randn(length(x),2));
   
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
        save(varargin{1},'W_eval','w_poly_fnc','dw_poly_x1_fnc','W_upper');
        
    end
end