figure(); 
subplot(2,1,1)
hold on

%RCI set
for i = 1:round(0.05*Tp/dt):length(MP_state)
    Ellipse_plot(M_ccm_ang,MP_state(i,[1,3])',25,'k',0.1);
end
%Terminal set
Ellipse_plot(P(1:2,1:2)*(1/alpha),x_eq([1;3]),25,'r');

%Full MP traj
plot(MP_state(:,1),MP_state(:,3),'r--','linewidth',1);
plot(test_state(1),test_state(3),'go','markersize',10,'markerfacecolor','g');
grid on
axis equal
xlabel('q'); ylabel('$\theta$','interpreter','latex');
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)


subplot(2,1,2)
hold on

%RCI set
for i = 1:round(0.05*Tp/dt):length(MP_state)
    Ellipse_plot(M_ccm_omg,MP_state(i,[2,4])',25,'k',0.1);
end
%Terminal set
Ellipse_plot(P(1:2,1:2)*(1/alpha),x_eq([2;4]),25,'r');

%Full MP traj
plot(MP_state(:,2),MP_state(:,4),'r--','linewidth',1);
plot(test_state(2),test_state(4),'go','markersize',10,'markerfacecolor','g');
grid on
axis equal
xlabel('$\dot{q}$','interpreter','latex'); ylabel('$\dot{\theta}$','interpreter','latex');
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)