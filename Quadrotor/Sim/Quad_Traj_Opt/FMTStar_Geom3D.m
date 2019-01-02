function [FMT_time, path_final] = FMTStar_Geom3D(N,V,N_nn,C_nn,EPS,obs,x0,goal,fig)

%N: num of nodes
%V: N x 4 node list
%N_nn: nearest neighbor index graph
%C_nn: nearest neighbor cost graph
%EPS: connection radius
%obs: collision obstacles
%x0: initial start location
%goal: goal polyhedron

%% Augment Nodes

%add start node to list
V = [x0;V];
N = N+1;

%augment graph matrices
N_nn = blkdiag(0,N_nn);
C_nn = blkdiag(0,C_nn);

%get connections from start node
dists = norms(V(2:N,:)-repmat(V(1,:),N-1,1),2,2);
N_nn(1,2:N) = (dists <= EPS)';
C_nn(1,2:N) = N_nn(1,2:N).*dists';

%% Initialize

%unvisited nodes (all except start)
V_unvisited = ones(N,1); 
V_unvisited(1) = 0;

%open nodes (none except start)
V_open = zeros(N,1); V_open(1) = 1;
C_open = Inf(N,1); C_open(1) = 0;

%graph nodes parent node indexes
V_graph = zeros(N,1);

%expand node (idx)
z = 1;
z_in_goal = 0;

%% Search

tic
while ~z_in_goal && sum(V_open) > 0
    
    %fwd reachable nodes of z
    NF_z = find(N_nn(z,:)); %idx of nn(z)
    
    for i = 1:length(NF_z)
        x = NF_z(i);
        
        %get unvisited neighbors of z
        if V_unvisited(x) 
            
            %bwd reachable nodes of x
            BN_x = find(N_nn(:,x)); %idx of nn(x)
            
            %find lowest cost node in V_open to x
            c_min = Inf; y_min = Inf;
            for j = 1:length(BN_x)
                y = BN_x(j);
                
                if V_open(y) 
                    c = C_open(y) + C_nn(y,x);
                    if c < c_min
                        c_min = c;
                        y_min = y; %store min node index
                    end
                end
            end
            if (y_min < Inf) && checkCollision(V(y_min,:),V(x,:),obs)
                %add x to graph
                V_graph(x) = y_min; 
                %add x to V_open
                V_open(x) = 1; C_open(x) = c_min;
                %remove x from V_unvisited
                V_unvisited(x) = 0;
            end
        end
    end

    %remove z from V_open
    V_open(z) = 0;
    C_open(z) = Inf;
    
    %find next min cost node in V_open
    [~,z] = min(C_open);
    
    %check if in goal
    z_in_goal = goal.contains(V(z,:)');
    
end
FMT_time = toc;

if (~z_in_goal)
    %exited with failure
    disp('Failure');
    path_final = 0;
    return;
end

%% Recover optimal path

% Search backwards from goal to start to find the optimal least cost path
path_idx = []; path_idx(1) = z; %start at end
q = z;
while (V_graph(q)~=0)
   %get parent
   p = V_graph(q);
   path_idx = [path_idx; p];
   q = p;   
end
%flip
path_idx = flipud(path_idx); %path node indices

%path waypoints
path = V(path_idx,:);

%Do plot
figure(fig)
hp1 = plot3(path(:,1),path(:,2),path(:,3),'bo-','linewidth',1.5);
drawnow
pause(1)


%% Smooth path

path = triangle_reduce(path,obs);

figure(fig)
hp2 = plot3(path(:,1),path(:,2),path(:,3),'ro-','linewidth',1.5);
drawnow
pause(1)
delete(hp1)

n_wp = size(path,1);
N_smooth_it = 1;
for i = 1:N_smooth_it
    
    %Cut corners
    path_up(1,:) = path(1,:);
    idx_a = 1; %count of updated path
    idx = 2; %corner node
    while (idx <= n_wp-1)
        shortcut = cut_corner(path_up(idx_a,:),path(idx,:),path(idx+1,:),...
                                     obs);
        path_up(idx_a+1:idx_a+2,:) = shortcut;
        idx_a = idx_a + 2;
        idx = idx+1;        
    end
    path_up(idx_a+1,:) = path(idx,:);
    
    %smooth
    path = triangle_reduce(path_up,obs);
    n_wp = size(path,1);
    
    %update plot
    figure(fig)
    delete(hp2)
    hp2 = plot3(path(:,1),path(:,2),path(:,3),'bo-','linewidth',1.5);
    drawnow;
    pause(1);
end

path_final = path;

end

function corner = cut_corner(p1,p2,p3,obstacles)
m1 = 0.5*(p1+p2);
m2 = 0.5*(p2+p3);
while ~checkCollision(m1,m2,obstacles)
    m1 = 0.5*(m1 + p2);
    m2 = 0.5*(m2 + p2);
end
corner = [m1;m2];
end
