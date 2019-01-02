
%% Create Tree obstacles (ll/ur format)

halton_pos = generateHaltonSamples(2,100);

%location of the trees
numTrees = 10;
tree_pos_xy = halton_pos(randperm(100,numTrees),:);
tree_pos_xy = tree_pos_xy.*repmat(World_dim(1:2),numTrees,1);

pine_rad = 1/4;
pine_height = 5/4;
tree_rad = 1.5/4;
tree_height = 7/4;
numTop = 6;
numDisc = 4;

tree_bounds = [];
tree_obstacles = [];

for i = 1:numTrees
    if rand(1) > 0.5
        tree_obstacles = [tree_obstacles;
                         createPineTreeObs([tree_pos_xy(i,:) 0], pine_rad, pine_height, numTop, numDisc)];
        tree_bounds = [tree_bounds;
                     tree_pos_xy(i,1)-pine_rad,tree_pos_xy(i,2)-pine_rad,0;
                     tree_pos_xy(i,1)+pine_rad,tree_pos_xy(i,2)+pine_rad,pine_height];
    else
        tree_obstacles = [tree_obstacles;
                         createTreeObs([tree_pos_xy(i,:) 0],tree_rad, tree_height, numTop, numDisc)];
        tree_bounds = [tree_bounds;
                     tree_pos_xy(i,1)-tree_rad,tree_pos_xy(i,2)-tree_rad,0;
                     tree_pos_xy(i,1)+tree_rad,tree_pos_xy(i,2)+tree_rad,tree_height];
    end
end

%% Create towers (ll/ur format)

numTowers = 6;
tower_pos_xy = halton_pos(randperm(100,numTowers),:);
tower_width = (9/4-0.1)*(World_dim(1)/15); tower_height = World_dim(3);
tower_pos_xy = tower_pos_xy.*repmat(World_dim(1:2)-tower_width,numTowers,1);

tower_obstacles = [];

for i = 1:numTowers
    tower = createBoxObs([tower_pos_xy(i,1),tower_pos_xy(i,2),0],[tower_width,tower_width,tower_height]);
    tower_obstacles = [tower_obstacles; tower];
end

%% Inflate

%tightest bounding box form
obstacles = [tree_bounds; tower_obstacles];
n_obstacles = size(obstacles,1)/2;


tube_infl = sqrt(diag(M_ccm_pos\eye(3)));
tube_infl = tube_infl';

size_infl = [0.2, 0.2, 0.2];

%add inflation by size of quad
obstacles_infl = obstacles  - kron(repmat([size_infl(1:2),0],n_obstacles,1),[1;0]) +...
                              kron(repmat(size_infl,n_obstacles,1),[0;1]);                       

%add inflation by size of tube bound
obstacles_coll = obstacles_infl  - kron(repmat([tube_infl(1:2),0],n_obstacles,1),[1;0]) +...
                                   kron(repmat(tube_infl,n_obstacles,1),[0;1]);
%% Plot 

fig = figure();
plot_all_obstacles();
plot3dObstacles(obstacles_coll,'r',0.0,1);

