%% Plot
close all

%% State Trajectory
figure()
hold on
plot(solve_t, X,'linewidth',2);
grid on
xlabel('Time [s]'); ylabel('States'); 
h_leg = legend('$q$','$\dot{q}$','$\theta$','$\dot{\theta}$');
set(h_leg,'interpreter','latex');
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)

%% Control Trajectory
figure()
hold on
plot(0:dt:t_end,MP_ctrl,'b-','linewidth',2);
plot(0:dt:t_end,U,'r-','linewidth',2);
xlabel('Time [s]');
ylabel('u(t)'); 
legend('nom','actual');
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)
grid on

figure()
plot(solve_t(1:end-1),U_fb,'linewidth',2); 
xlabel('Time [s]');
ylabel('k(x^{*},x)');
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)
grid on

%% 2D State Plot
visualize_FLR;
subplot(2,1,1);
plot(X(:,1),X(:,3),'b-','linewidth',2);
subplot(2,1,2);
plot(X(:,2),X(:,4),'b-','linewidth',2);

%% Solve Time
figure()
plot(solve_t(1:end-1),ctrl_solve_time,'rd','markersize',10,'markerfacecolor','k');
grid on
legend('CCM');
xlabel('Time [s]');
ylabel('Solve time [s]'); title('Solve time');
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)

%% Solve success
figure()
plot(solve_t(1:end-1),opt_solved,'rd','markersize',10,'markerfacecolor','k');
grid on
legend('CCM (0,1,6)');
xlabel('Time [s]');
title('Convergence');
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)
 
%% Geod distances

figure()
plot(solve_t(1:end-1),sqrt(geo_energy),'b-','linewidth',2);
hold on
plot(solve_t(1:end-1),d_bar*ones(T_steps,1),'r-','linewidth',2);
grid on
legend('d(x^{*},x)','RCI bound');
xlabel('Time [s]');
    set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)