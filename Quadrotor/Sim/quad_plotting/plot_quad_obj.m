function [h_quad,R] = plot_quad_obj(r,rot,varargin)

roll = rot(1);
pitch = rot(2);
yaw = rot(3);

R_z = [cos(yaw), -sin(yaw), 0;
       sin(yaw), cos(yaw), 0;
       0,0,1];
   
R_y = [cos(pitch), 0, sin(pitch);
       0, 1, 0;
      -sin(pitch), 0, cos(pitch)];
  
R_x = [1,0,0;
       0, cos(roll), -sin(roll);
       0, sin(roll), cos(roll)];
  
R = R_x*R_y*R_z;

%% distort and translate

%quadrotor body 
bound = [0.1,0.1, 0;
           0.1, -0.1, 0;
           -0.1,-0.1,0;
           -0.1,0.1,0]';

%rotor base       
N_circ = 20;       
base_circle = circle([0;0;0],0.05,N_circ); 

%axes base
ax = 0.3*eye(3);

%NED -> NWU
Flip_M = [1,0,0;
          0,-1,0;
          0,0,-1];
      
%define points
q_bound = bound;
q_rotors = cell(4,1);
q_ax = Flip_M*R*ax + kron(ones(1,3),Flip_M*r);

c_bound = Flip_M*R*base_circle;
for j = 1:4
    q_bound(:,j) = Flip_M*r + Flip_M*R*bound(:,j);
    
    q_rotors{j} = kron(ones(1,N_circ),q_bound(:,j)) + ...
                   c_bound;
end

%% plot

if nargin > 2
    figure(varargin{1})
else
    figure()
end


%draw objects
h_rotors = zeros(4,1);
for j = 1:4
    h_rotors(j) = patch(q_rotors{j}(1,:),q_rotors{j}(2,:),q_rotors{j}(3,:),'k','FaceAlpha',0.1,'linewidth',2);
end

h_ax = zeros(3,1);
col_ax = {'r','b','g'};
for j = 1:3
    h_ax(j) = line([r(1);q_ax(1,j)],[-r(2),q_ax(2,j)],[-r(3),q_ax(3,j)],...
                     'color',col_ax{j},'linewidth',2);
end

h_arms = zeros(2,1);
for j = 1:2
    h_arms(j) = line([q_bound(1,j);q_bound(1,j+2)],...
                     [q_bound(2,j);q_bound(2,j+2)],...
                     [q_bound(3,j);q_bound(3,j+2)],...
                     'color','k','linewidth',2);
end

axis equal

%% return

h_quad = [h_rotors;h_ax;h_arms];

end