
%% Trajectory plot (NWU)
figure();
plot_all_obstacles();
plot3(X(:,1),-X(:,2),-X(:,3),'b-','linewidth',2); hold on
plot3(MP_state(:,1),-MP_state(:,2),-MP_state(:,3),'r-','linewidth',2);
plot3(X(1,1),-X(1,2),-X(1,3),'go','markersize',15,'markerfacecolor','g');
grid on; axis tight;
xlabel('x'); ylabel('y'); zlabel('h');
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)

%% Trajectory errors
figure()
subplot(2,1,1)
plot(solve_t(1:end-1),err(:,1:3),'linewidth',2);
grid on
xlabel('Time [s]');
ylabel('[m]');
legend('e_x','e_y','e_z');
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)

subplot(2,1,2)
plot(solve_t(1:end-1),err(:,4)*(180/pi),'linewidth',2);
grid on
xlabel('Time [s]'); ylabel('[deg]');
legend('e_\psi');
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)

%% Velocity
figure()
plot(MP_t(1:X_steps),X(1:end,4:6),'linewidth',2);
grid on
xlabel('Time [s]');
ylabel('Velocity [m/s]');
legend('x','y','z');
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)

%% Control effort
figure()
plot(MP_t, MP_ctrl(:,1),'--','linewidth',2); hold on
plot(MP_t(1:X_steps) ,U(:,1),'-','linewidth',2);
xlabel('Time [s]');
legend('nom','actual');
ylabel('$\dot{\tau}$ [1/s]','interpreter','latex');
grid on

figure()
subplot(3,1,1)
plot(MP_t, MP_ctrl(:,2),'--','linewidth',2); hold on
plot(MP_t(1:X_steps) ,U(:,2),'-','linewidth',2);
xlabel('Time [s]');
legend('nom','actual');
ylabel('$\dot{\phi}$ [rad/s]','interpreter','latex');
grid on
% set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)

subplot(3,1,2)
plot(MP_t, MP_ctrl(:,3),'--','linewidth',2); hold on
plot(MP_t(1:X_steps), U(:,3),'-','linewidth',2);
xlabel('Time [s]');
legend('nom','actual');
ylabel('$\dot{\theta}$ [rad/s]','interpreter','latex');
grid on
% set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)

subplot(3,1,3)
plot(MP_t, MP_ctrl(:,4),'--','linewidth',2); hold on
plot(MP_t(1:X_steps), U(:,4),'-','linewidth',2);
xlabel('Time [s]');
legend('nom','actual');
ylabel('$\dot{\psi}$ [rad/s]','interpreter','latex');
grid on
% set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)

%% Rotation angles
figure()
subplot(2,1,1)
plot(MP_t, MP_state(:,8)*(180/pi),'--','linewidth',2); hold on
plot(MP_t(1:X_steps), X(:,7)*(180/pi),'-','linewidth',2); hold on
xlabel('Time [s]');
ylabel('$\phi$ [deg]','interpreter','latex');
grid on
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)

subplot(2,1,2)
plot(MP_t, MP_state(:,9)*(180/pi),'--','linewidth',2); hold on
plot(MP_t(1:X_steps), X(:,8)*(180/pi),'-','linewidth',2); hold on
xlabel('Time [s]');
ylabel('$\theta$ [deg]','interpreter','latex');
grid on
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)

%% Geodesic Energy
figure()
plot(solve_t(1:end-1),geo_energy,'b-','linewidth',2); hold on
plot(solve_t(1:end-1),(d_bar^2)*ones(length(solve_t)-1,1),'k-','linewidth',2);
grid on
xlabel('Time [s]'); ylabel('Energy');
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)

%% Computation time

figure()
plot(solve_t(1:end-1),ctrl_solve_time,'rd','markersize',10,'markerfacecolor','k');
grid on
xlabel('Time [s]');
ylabel('Solve time [s]'); title('Solve time');
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)

%% Solve success

figure()
hold on
plot(solve_t(1:end-1),opt_solved,'rd','markersize',10,'markerfacecolor','k');
grid on
xlabel('Time [s]');
title('Convergence');
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)
 
%% Animate
keyboard;
% close all;
plot_quad_movie(MP_state(:,1),MP_state(:,2),MP_state(:,3),MP_t,X,round(0.1/dt_sim),M_ccm_pos,goal,...
                tree_obstacles,tower_obstacles,obstacles_infl,World_dim);
