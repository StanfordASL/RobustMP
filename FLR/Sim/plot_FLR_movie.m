function plot_FLR_movie(t,x_act,t_step)
%t: time history
%x_act: state history
%t_step: time-stepping (ints w.r.t. length(t))

motor_bound = [[-0.15;-0.05],[-0.1;0.5],[0.1;0.5],[0.15;-0.05]];
link_bound = [[-0.1;0.0],[-0.1;1.4],[0.1;1.4],[0.1;0.0]];

link_point = [0;0.4];

figure();
plot(0,0,'ro','markersize',15,'markerfacecolor','r');
hold on;

xl = [-2, 2];
yl = [-2, 2];

i = 1;
motor_b = motor_bound;
link_b = link_bound;
l_p = link_point;
l_vec = 4*l_p;

R_mot = [cos(x_act(i,3)), sin(x_act(i,3));
        -sin(x_act(i,3)), cos(x_act(i,3))];
R_link = [cos(x_act(i,1)), sin(x_act(i,1));
        -sin(x_act(i,1)), cos(x_act(i,1))];
    
l_p = R_mot*l_p;
l_vec = R_link*l_vec + l_p;

for j = 1:4
    motor_b(:,j) = [0;0] + R_mot*motor_bound(:,j);
    link_b(:,j) = l_p + R_link*link_bound(:,j);
end
hm = patch(motor_b(1,:),motor_b(2,:),'r','FaceAlpha',0.1,'linewidth',2);
hm_l = line([0;2*l_p(1)],[0;2*l_p(2)],'color','r','linewidth',2);
hl = patch(link_b(1,:),link_b(2,:),'b','FaceAlpha',0.1,'linewidth',2);
hl_l = line([l_p(1);l_vec(1)],[l_p(2); l_vec(2)],'color','b','linewidth',2);
hold off

set(gca,'Xlim',xl);
set(gca,'Ylim',yl);

grid on;

axis manual;
keyboard;

for i = t_step:t_step:length(t)
    
    motor_b = motor_bound;
    link_b = link_bound;
    l_p = link_point;
    l_vec = 4*l_p;
    
    R_mot = [cos(x_act(i,3)), sin(x_act(i,3));
        -sin(x_act(i,3)), cos(x_act(i,3))];
    R_link = [cos(x_act(i,1)), sin(x_act(i,1));
        -sin(x_act(i,1)), cos(x_act(i,1))];
    
    l_p = R_mot*l_p;
    l_vec = R_link*l_vec + l_p;
    
    for j = 1:4
        motor_b(:,j) = [0;0] + R_mot*motor_bound(:,j);
        link_b(:,j) = l_p + R_link*link_bound(:,j);
    end
    
    hold on;
    set(hm,'XData',motor_b(1,:),'YData',motor_b(2,:));
    set(hm_l,'XData',[0;2*l_p(1)],'YData',[0;2*l_p(2)]);
    set(hl,'XData',link_b(1,:),'YData',link_b(2,:));
    set(hl_l,'XData',[l_p(1);l_vec(1)],'YData',[l_p(2);l_vec(2)]);
    hold off;
    
    title(sprintf('t = %.2f',t(i)));
    drawnow;
       
    
end
end