function dx_dot = ode_sim(t,x,t_span,u_nom,k,f,B,B_w,w)

u = interp1(t_span,u_nom,t) + k;
dx_dot = f(x) + B(x)*u' + B_w*w;

end