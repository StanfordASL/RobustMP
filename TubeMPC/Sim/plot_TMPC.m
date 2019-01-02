%% Plot

close all;

%% 2D State Plot
figure()
for i_mpc = 1:T_steps_MPC
    plot(MPC_state{i_mpc}(1:round(delta/dt)+1,1),...
         MPC_state{i_mpc}(1:round(delta/dt)+1,2),'r-','linewidth',2);
      hold on;
    Ellipse_plot(M_ccm*(1/d_bar^2),MPC_state{i_mpc}(1,1:2),25,'g',0.5);
    Ellipse_plot(M_ccm*(1/d_bar^2),MPC_state{i_mpc}(1+round(delta/dt),1:2),25,'r',0.5);
end
plot(X(:,1),X(:,2),'b-','linewidth',2);    

Ellipse_plot((1/alpha)*P,x_eq,25,'r');
xlabel('x_1'); ylabel('x_2');
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)
grid on; 

%% State Trajectory
figure()
hold on
plot(solve_t, X,'linewidth',2);
grid on
xlabel('Time [s]'); ylabel('States'); 
h_leg = legend('$x_1$','$x_2$');
set(h_leg,'interpreter','latex');
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)

%% Control Trajectory
figure()
hold on
plot(0:dt:t_end,U_nom,'b-','linewidth',2);
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

%% Solve Time
figure()
hold on
plot(solve_t(1:end-1),ctrl_solve_time(:,1),'ro','markersize',10,'markerfacecolor','g');
plot(solve_t(1:end-1),ctrl_solve_time(:,2),'rd','markersize',10,'markerfacecolor','k');
grid on
legend('MPC','CCM');
xlabel('Time [s]');
ylabel('Solve time [s]'); title('Solve time');
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)

%% Solve success
figure()
hold on
plot(solve_t(1:end-1),opt_solved(:,1),'ro','markersize',15,'markerfacecolor','g');
plot(solve_t(1:end-1),opt_solved(:,2),'rd','markersize',10,'markerfacecolor','k');
grid on
legend('MPC (1,2,3)','CCM (0,1,6)');
xlabel('Time [s]');
title('Convergence');
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)
 
%% Geod distances
figure()
plot(solve_t(1:end-1),sqrt(geo_energy(:,1)),'b-','linewidth',2);
hold on
plot(solve_t(1:end-1),sqrt(geo_energy(:,2)),'bo','markersize',10,'markerfacecolor','b','linewidth',2);
plot(solve_t(1:end-1),d_bar*ones(length(solve_t)-1,1),'r-','linewidth',2);
grid on
legend('d(x^{*},x)','RCI bound');
xlabel('Time [s]');
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)