%2D State Plot
% P_inv = diag([39.0251,486.0402]);
figure()
plot(MP_state(:,1),MP_state(:,2),'r-','linewidth',2); hold on
Ellipse_plot((1/alpha)*P,x_eq,25,'r');
xlabel('x_1'); ylabel('x_2');
set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)
grid on; 