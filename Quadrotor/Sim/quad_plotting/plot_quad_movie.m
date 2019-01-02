function plot_quad_movie(x_des,y_des,z_des,solve_t,X,~,M,goal,...
                         tree_obstacles,tower_obstacles,obstacles_infl,World_dim)

close all;
fig = figure();
plot3(x_des,-y_des,-z_des,'r-','linewidth',2); hold on
plot3(X(:,1),-X(:,2),-X(:,3),'k-','linewidth',2.0);
goal.plot('color','blue','alpha',0.3); 

%% FPV view:
fpv_view = 1;
set(gca,'CameraViewAngleMode','manual')
camproj('perspective')
camva(45)

%%
dt = solve_t(2)-solve_t(1);
t_step = round((1/20)/dt);
t_step_b = round(0.3/dt);
T_b = 2; %T_b second lookahead for bound
n_b = 1+round((T_b/0.3));

i = 1;
pos = X(i,1:3)';
att = X(i,7:9);

%% create objects

hold on
% if (~fpv_view)
    h_b = [];
    for j = 1:n_b
        j_idx = min(i + (j-1)*t_step_b, length(solve_t));
        h_b(j) = proj_Ellipse([1:3],M,1,[x_des(j_idx);-y_des(j_idx);-z_des(j_idx)],30,'b',0);
    end
% end
[h_quad,~] = plot_quad_obj(pos,att,fig);
rot = rot_matrix(att(1),att(2),atan2(X(i,5),X(i,4)));

plot3dObstacles(tree_obstacles,'g',0.6); hold on
plot3dObstacles(tower_obstacles,[0.5,0.5,0.5],1.0);
plot3dObstacles(obstacles_infl,'r',0.01,1);
plot3dObstacles([0 0 0; World_dim],'k',0);
patch([0,World_dim(1),World_dim(1),0],...
      [0, 0, World_dim(2), World_dim(2)],...
      [0, 0,  0, 0],'FaceColor',[0,0.5,0],'FaceAlpha',0.5);
  
hold off

view(3);

Flip_M = [1,0,0;0,-1,0;0,0,-1];

%set view
if (fpv_view)
    pos_tail = [-1;1;-1];
    pos_view = [0.1;0;0.0];
    c_pos = (Flip_M*pos)'+(Flip_M*rot*pos_tail)';
    campos(c_pos);
    c_target =(Flip_M*pos)'+(Flip_M*rot*pos_view)'; 
    camtarget(c_target);
else
    %iso-view
    campos((Flip_M*pos)'+[-1.5,-2.5,1]);
    camtarget((Flip_M*pos)');
end
drawnow;

%% Setup Movie
record_vid = 1;

if (record_vid)
    writerObj = VideoWriter('Quad_sim.mp4');
    writerObj.FrameRate = 1/(t_step*dt);
    writerObj.Quality = 100;
    open(writerObj);
    set(gcf, 'renderer', 'zbuffer')
end

%% Record

keyboard;
for i = 1:t_step:length(solve_t)
    pos = X(i,1:3)';
    att = X(i,7:9);
    
    delete(h_quad); 
    hold on
    %update quad and bound
%     if (~fpv_view)
        delete(h_b);
        h_b = [];
        for j = 1:n_b
            j_idx = min(i + (j-1)*t_step_b, length(solve_t));
            h_b(j) = proj_Ellipse([1:3],M,1,[x_des(j_idx);-y_des(j_idx);-z_des(j_idx)],30,'b',0);
        end
%     end
    [h_quad,~] = plot_quad_obj(pos,att,fig);
    rot = rot_matrix(att(1),att(2),atan2(X(i,5),X(i,4)));
    
    hold off

    %update view
    if (fpv_view)
        c_pos = (Flip_M*pos)'+(Flip_M*rot*pos_tail)';
        campos(c_pos);
        c_target =(Flip_M*pos)'+(Flip_M*rot*pos_view)';
        camtarget(c_target);
    else
        %iso-view
        campos((Flip_M*pos)'+[-1.5,-2.5,1]);
        camtarget((Flip_M*pos)');
    end
    drawnow;
    
    %record
    if (record_vid)
        thisFrame = getframe(gcf);
        %Write this frame out to a new video file.
        writeVideo(writerObj, thisFrame);
    end
end 

if (record_vid)
    close(writerObj);
end
end