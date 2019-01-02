%% 2D State Plot

close all; 

figure(1)
hold on
%Plot obstacles (inflated by vehicle size)
for i_ob = 1:obs.n_obs
    Ellipse_plot(eye(2)*(1/(obs_rad(i_ob)+len)^2),obs.pos(:,i_ob), 25,'k',1.0);
end

ellipse_color = [188,216,240]/255;
ellipse_alpha = 0.15;
plot_quivers = 0;

if (do_mpc)
    %if used time-varying bound for coll., then do time-varying Ellipsoid plot
    plot_time_var = obs_mpc.time_var;
    for i_mpc = 1:T_steps_MPC
        %MPC reference trajectory executed segement
        
        num_pnts = round(delta/0.05);
        
        %Outer geodesic ball at start of MPC segment
        E_start = geo_energy(1+(i_mpc-1)*(delta/dt_sim),2);
        if plot_time_var
            Ellipse_plot(M_ccm_pos_unscaled*(1/E_start),MPC_state{i_mpc}(1,1:2)',30,'g',ellipse_alpha);
        else
            Ellipse_plot(M_ccm_pos,MPC_state{i_mpc}(1,1:2)',30,ellipse_color,ellipse_alpha);
        end
        
        t_mpc_span = linspace(0,delta,num_pnts);
        %Evolution of outer geodesic ball over MPC segment
        for j = 2:length(t_mpc_span)-1
            E_j = (sqrt(E_start)*exp(-lambda*t_mpc_span(j)) + d_bar*(1-exp(-lambda*t_mpc_span(j))))^2;
            if plot_time_var
                Ellipse_plot(M_ccm_pos_unscaled*(1/E_j),MPC_state{i_mpc}(round(t_mpc_span(j)/dt)+1,1:2)',30,[0,1,0]+(j/length(t_mpc_span))*[0,-1,1],ellipse_alpha);
            else
                Ellipse_plot(M_ccm_pos,MPC_state{i_mpc}(round(t_mpc_span(j)/dt)+1,1:2)',30,ellipse_color,ellipse_alpha);
            end
        end
        %Final outer geodesic ball
        E_end = (sqrt(E_start)*exp(-lambda*delta) + d_bar*(1-exp(-lambda*delta)))^2;
        if plot_time_var
            Ellipse_plot(M_ccm_pos_unscaled*(1/E_end),MPC_state{i_mpc}(round(delta/dt)+1,1:2)',30,'r',ellipse_alpha);
        else
            Ellipse_plot(M_ccm_pos,MPC_state{i_mpc}(round(delta/dt)+1,1:2)',30,ellipse_color,ellipse_alpha);
        end
        
        if i_mpc == T_steps_MPC
            plot(MPC_state{i_mpc}(:,1),MPC_state{i_mpc}(:,2),'r--','linewidth',1.5);
        else
            plot(MPC_state{i_mpc}(1:round(delta/dt)+1,1),MPC_state{i_mpc}(1:round(delta/dt)+1,2),'r--','linewidth',1.5);
        end
        
    end
else
    for i = 1:(0.5/dt):length(MP_t)
        Ellipse_plot(M_ccm_pos,MP_state(i,1:2),30,'k',0.15);
    end
end

if (plot_quivers)
    quiver(X(1:(0.1/dt):end,1),X(1:(0.1/dt):end,2),...
        -sin(X(1:(0.1/dt):end,3)),...
         cos(X(1:(0.1/dt):end,3)),0.5);
end

%Nominal motion plan
plot(MP_state(:,1),MP_state(:,2),'k--','linewidth',1.5);

%Plot actual trajectory
plot(X(:,1),X(:,2),'m-','linewidth',2);

%Final condition
Ellipse_plot(P(1:2,1:2)*(1/(alpha)), x_eq(1:2),30,'k');
xlabel('$X$','interpreter','latex'); 
ylabel('$Z$','interpreter','latex');
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)
grid off; 
% axis equal

%% State Trajectory
figure()
hold on
plot(solve_t, X(:,3:6),'linewidth',2);
grid on
xlabel('Time [s]'); ylabel('States'); 
h_leg = legend('$\phi$','$v_x$','$v_z$','$\dot{\phi}$');
set(h_leg,'interpreter','latex');
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)

%% Control Trajectory
figure()
subplot(2,1,1)
hold on
plot(0:dt:t_end,U_nom(:,1),'b-','linewidth',2);
plot(0:dt:t_end,U(:,1),'r-','linewidth',2);
xlabel('Time [s]');
ylabel('u(t)'); 
legend('nom','actual');
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)
grid on

subplot(2,1,2)
hold on
plot(0:dt:t_end,U_nom(:,2),'b-','linewidth',2);
plot(0:dt:t_end,U(:,2),'r-','linewidth',2);
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
if (do_mpc)
%     plot(solve_t(1:end-1),d_bar*ones(T_steps,1),'r-','linewidth',2);
    for i_mpc = 1:T_steps_MPC
        E_start = geo_energy(1+(i_mpc-1)*(delta/dt_sim),2);
        t_start = solve_t(1+(i_mpc-1)*(delta/dt_sim));
        t_mpc_span = 0:dt_sim:delta;
        d_mpc = sqrt(E_start)*exp(-lambda*t_mpc_span) + d_bar*(1-exp(-lambda*t_mpc_span));
        plot(t_mpc_span+t_start,d_mpc,'r-','linewidth',2);
    end
    plot(solve_t(1:end-1),sqrt(geo_energy(:,2)),'bo','markersize',10,'markerfacecolor','b','linewidth',2);
else
    plot(solve_t(1:end-1),sqrt(geo_energy(1,1))*exp(-lambda*solve_t(1:end-1)) + d_bar*(1-exp(-lambda*solve_t(1:end-1))),'r-','linewidth',2);
end
grid on
legend('d(x^{*},x)','RCI bound');
xlabel('Time [s]');
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)

%% MPC Rejoin time

if (do_mpc)
    figure()
    plot(MPC_time(:,1),MPC_time(:,2),'ro','markersize',10,'markerfacecolor','r');
    hold on
    plot(MPC_time(:,1),min([MPC_time(:,1)+delta,Tp*ones(T_steps_MPC,1)],[],2),'bo','markersize',10,'markerfacecolor','b');
    grid on
    xlabel('Time [s]');ylabel('t_i + T');
    legend('Optimal re-join time','Minimum re-join time');
    set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)
end