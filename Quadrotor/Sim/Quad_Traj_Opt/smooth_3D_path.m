function [pos_coef,vel_coef,acc_coef,jer_coef,T_seg,coeff] = smooth_3D_path(path,obs,fig)
%path: FMT* (or waypoint) path

%smoothing optimization:
%J = p^T Q p + k_T*sum(T_i)
%s.t. A*p = d = [d_F;d_P]; where d_F: known, d_P: unknown
%set p = A\d <-> J = d^T (A^-T Q A^-1) d + sum(T_i)

%Fixed T_i, optimize for d: get_dP function (1 step exact)
%Fixed d, optimize for T: get_T function (n_T steps gradient descent)
%alternate

v_init = zeros(1,3);

k_T = 100;
n_T = 1; %n_T gradient steps for T optimization

%% Setup polynomial structure

%order
p_order = 10;

%coefficient return for pos,vel,acc,jerk
pos_coef = @(T) [1, T, T^2, T^3, T^4, T^5, T^6, T^7, T^8, T^9];
vel_coef = @(T) [0, 1, 2*T, 3*T^2, 4*T^3, 5*T^4, 6*T^5, 7*T^6, 8*T^7, 9*T^8];
acc_coef = @(T) [0, 0, 2, 6*T, 12*T^2, 20*T^3, 30*T^4, 42*T^5, 56*T^6, 72*T^7];
jer_coef = @(T) [0, 0, 0, 6, 24*T, 60*T^2, 120*T^3, 210*T^4,336*T^5,504*T^6];
snap_coef = @(T) [0,0,0,0,24,120*T,360*T^2,840*T^3,1680*T^4,3024*T^5];

%% Get Min Snap Q matrix for single T segment

p = sym('p',[p_order,1]);
syms T
snap_poly_sq = (snap_coef(T)*p).^2;
J = int(snap_poly_sq,T);

Q = sym('Q',[p_order,p_order]);
Q(1:4,:) = zeros(4,p_order);
for i = 5:p_order
    for j = i:p_order
        Q(i,j) = diff(diff(J,p(i)),p(j))/2;
    end
end
for i = 5:p_order
    for j = 1:i-1
        Q(i,j) = Q(j,i);
    end
end
Q_fun = matlabFunction(Q); %symbolic function

%% Start

done = 0;

while (~done)
    N_nodes = size(path,1);
    N_seg = N_nodes-1;
    
    %% Get Full Gram matrices
    
    Tv = sym('Tv',[N_seg,1]);
    Q_full = Q_fun(Tv(1));
    for ns = 2:N_seg
        Q_full = blkdiag(Q_full,Q_fun(Tv(ns)));
    end
    
    Q_full_fun = matlabFunction(Q_full,'Vars',{Tv});
    
    A = get_constraints(N_seg,Tv,p_order,pos_coef,vel_coef,acc_coef,jer_coef,snap_coef);
    A_inv = inv(A);
%     G = A_inv.'*Q_full*A_inv;
%     %now convert to matlab functions
%     G_T_full = cell(N_seg,1);
%     G_T_fun = cell(N_seg,1);
%     for ns = 1:N_seg
%         G_T_full{ns} = diff(G,Tv(ns));
%         G_T_fun{ns} = matlabFunction(G_T_full{ns},'Vars',{Tv});
%     end
    
    A_inv_fun = matlabFunction(A_inv,'Vars',{Tv});
    
    %% Setup d_F
    
    v0 = v_init;
    
    d_wp = zeros(2*N_seg,3);
    for ns = 1:N_seg
        d_wp(1+(ns-1)*2,:) = path(ns,:);
        d_wp(2*ns,:) = path(ns+1,:);
    end
    d_F = [d_wp;
        v0;
        zeros(2+2+4*(N_seg-1),3)];
    s_F = size(d_F,1);
    
    %% Initial guess for T_seg
    
    T_seg = zeros(N_seg,1);
    %guess some segment times 
    for ns = 1:N_seg
        T_seg(ns) = norm(path(ns+1,:)-path(ns,:))/3.0;
    end
    
    %Do min snap
    A_inv = A_inv_fun(T_seg);
    d_P = get_dP(d_F,s_F,T_seg,Q_full_fun,A_inv);
    
    %Get J
    d = [d_F;d_P];
    J = trace(d'*A_inv'*Q_full_fun(T_seg)*A_inv*d) + k_T*sum(T_seg);
    
    %% Begin alternating algorithm
    
    max_it = 10000;
    converged = 1;
    eps = 1e-3;
    step_T = 1e-3;
    itr = 1;
    J_vec = J*ones(max_it,1);
    
    while (~converged) && (itr<max_it)
        
        %% Min T
        T_seg_up = get_T(G_T_fun,k_T,T_seg,d_F,d_P,n_T,step_T/itr);
        
        %% Min snap
        A_inv = A_inv_fun(T_seg_up);
        d_P = get_dP(d_F,s_F,T_seg_up,Q_full_fun,A_inv);
        
        %% Evaluate cost
        d = [d_F;d_P];
        J_up = trace(d'*A_inv'*Q_full_fun(T_seg_up)*A_inv*d) + k_T*sum(T_seg_up);
        
        %% Check convergence
        dJ = (J_up-J)/J;
        if abs(dJ)<=eps
            converged = 1;
        end
        
        %% Update
        disp([T_seg,T_seg_up]);
        fprintf('%d, %.4f, %.4f, %.4f \n',itr,J,J_up,dJ);
        itr = itr + 1;
        J_vec(itr) = J_up;
        
        J = J_up;
        T_seg = T_seg_up;
        
    end
    
%     figure()
%     plot(J_vec(1:itr),'linewidth',2);
%     ylabel('Smoothing cost');
%     grid on
        
    %% Extract optimal solution
    
    A_inv = A_inv_fun(T_seg);
    d_P = get_dP(d_F,s_F,T_seg,Q_full_fun,A_inv);
    coeff = A_inv*[d_F;d_P];
    
    %% Compute position polynomials
    
    poly_seg = cell(N_seg,1);
    dt = 0.1;
    for ns = 1:N_seg
        len = floor(T_seg(ns)/dt)+1;
        poly_seg{ns} = zeros(len,3);
        for i = 1:len
            t = (i-1)*dt;
            poly_seg{ns}(i,:) = pos_coef(t)*coeff(1+(ns-1)*p_order:ns*p_order,:);
        end
    end
    
    %% Plot
    
    figure(fig)
    hp = [];
    for ns = 1:N_seg
        hp(ns) = plot3(poly_seg{ns}(:,1),poly_seg{ns}(:,2),poly_seg{ns}(:,3),'b-','linewidth',2);
    end
    pause(1);
    
    %% Collision check polyspline
    
    disp('Checking for collisions...');
    detected_coll = 0;
    for ns = 1:N_seg
        fprintf('seg:%d\n',ns);
        %check segment ns
        for i = 1:size(poly_seg{ns},1)-1
            if (detected_coll)
                disp('found');
                break;
            end
            detected_coll = ~checkCollision(poly_seg{ns}(i,:),...
                poly_seg{ns}(i+1,:),...
                obs);
        end
        %if found collision in segment ns
        if (detected_coll)
            %add waypoint between nodes ns and ns+1
            wp_new = 0.5*(path(ns,:)+path(ns+1,:));
            path_up = [path(1:ns,:);
                wp_new;
                path(ns+1:end,:)];
        end
        if (detected_coll)
            break;
        end
    end
    
    if (~detected_coll)
        disp('no collision found');
        break;
    else
        disp('collision found');
        path = path_up;
        delete(hp);
        plot3(path(:,1),path(:,2),path(:,3),'ro-','linewidth',1.5);
    end    
end

end



