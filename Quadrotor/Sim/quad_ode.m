function dx_dot = quad_ode(t,x,t_span,u_nom,k,thrust,f,B,B_w,w)

global kp_om Jq;

uc = interp1(t_span,u_nom,t) + k; %roll_dot,pitch_dot,yaw_dot

%get thrust
% thrust = x(end);

%get torque
euler_dot_des = uc';
om_des = R_om(x(7:9))*euler_dot_des;
om = x(10:12);
M = kp_om*(om_des - om);

%assemble
u = [thrust;M];

%propagate
dx_dot = f(x) + B(x)*u + B_w*w;

end