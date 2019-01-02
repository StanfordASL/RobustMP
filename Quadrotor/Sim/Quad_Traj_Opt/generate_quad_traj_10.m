function [t_vec,xc_nom,uc_nom,thrust_nom] = generate_quad_traj_10(dt,mq,g,pos_init,varargin)

sim_case = 'polyspline';

switch (sim_case)
    case 'fig8'
        T = 10;
        t_vec = (0:dt:T)';
        om = 2*pi*0.125;
        x_d = 2*cos(om*t_vec);
        y_d = 3*0.5*sin(2*om*tt_vec);
        z_d = zeros(length(t_vec),1);
        yaw_d = zeros(length(t_vec),1);
        
        vx_d = -2*om*sin(om*t_vec);
        vy_d =  3*om*cos(2*om*t_vec);
        vz_d = zeros(length(t_vec),1);
        yd_d = zeros(length(t_vec,1));
        
        ax_d = -2*(om^2)*cos(om*t_vec);
        ay_d = -3*2*(om^2)*sin(2*om*t_vec);
        az_d =  zeros(length(t_vec),1);
        
        jx_d =  2*(om^3)*sin(om*t_vec);
        jy_d = -3*4*(om^3)*cos(2*om*t_vec);
        jz_d =  zeros(length(t_vec),1);
        
    case 'circle'
        T = 20;
        t_vec = (0:dt:T)';
        radius = 1;
        
        om = 2*3.14*(1.0/10);
        x_d = pos_init(1)+radius*sin(om*t_vec);
        z_d = pos_init(3)*ones(length(t_vec),1);
        y_d = pos_init(2) + radius-radius*cos(om*t_vec);
        yaw_d = zeros(length(t_vec),1);
        
        vx_d = radius*om*cos(om*t_vec);
        vz_d = zeros(length(t_vec),1);
        vy_d = radius*om*sin(om*t_vec);
        yd_d = zeros(length(t_vec),1);
        
        ax_d = -radius*(om^2)*sin(om*t_vec);
        az_d = zeros(length(t_vec),1);
        ay_d = radius*(om^2)*cos(om*t_vec);
        
        jx_d = -radius*(om^3)*cos(om*t_vec);
        jz_d = zeros(length(t_vec),1);
        jy_d = -radius*(om^3)*sin(om*t_vec);
        
    case 'helix'
        T = 10;
        t_vec = (0:dt:T)';
        
        om = 2*pi*0.125;
        x_d = sin(om*t_vec);
        y_d = cos(om*t_vec);
        z_d = -0.1*t_vec;
        yaw_d = zeros(length(t_vec),1);
        
        vx_d = om*cos(om*t_vec);
        vy_d = -om*sin(om*t_vec);
        vz_d = -0.1*ones(length(t_vec),1);
        yd_d = zeros(length(t_vec),1);
        
        ax_d = -(om^2)*sin(om*t_vec);
        ay_d = -(om^2)*cos(om*t_vec);
        az_d = zeros(length(t_vec),1);
        
        jx_d = -(om^3)*cos(om*t_vec);
        jy_d =  (om^3)*sin(om*t_vec);
        jz_d = zeros(length(t_vec),1);
        
    case 'hover'
        T = 5;
        tau = T/4;
        t_vec = (0:dt:T)';
        
        x_d = pos_init(1)*ones(length(t_vec),1);
        y_d = pos_init(2)*ones(length(t_vec),1);
        z_d =  -1*(1 - 1*exp(-t_vec/tau)) + pos_init(3);
        yaw_d = zeros(length(t_vec),1);
        
        vx_d = zeros(length(t_vec),1);
        vy_d = zeros(length(t_vec),1);
        vz_d = -1*((1/tau)*exp(-t_vec/tau));
        yd_d = zeros(length(t_vec),1);
        
        ax_d = zeros(length(t_vec),1);
        ay_d = zeros(length(t_vec),1);
        az_d = -1*(-(1/tau^2)*exp(-t_vec/tau));
        
        jx_d = zeros(length(t_vec),1);
        jy_d = zeros(length(t_vec),1);
        jz_d = -1*((1/tau^3)*exp(-t_vec/tau));
        
    case 'polyspline'
        poly_file = varargin{1};
        load(poly_file);
        p_order = 10;
        
        %Create splines for each segment (also convert NWU to NED)
        roll_MAX = 50*(pi/180);
        roll_max = 60*(pi/180);
        slow_scale = 1.0;
        while (roll_max > roll_MAX)
            
            break_T = [0;cumsum(T_seg)];
            t_vec = (0:dt:break_T(end)*slow_scale)';
            x_d = zeros(length(t_vec),1); y_d = x_d; z_d = x_d;
            vx_d = x_d; vy_d = x_d; vz_d = x_d;
            ax_d = x_d; ay_d = x_d; az_d = x_d;
            jx_d = x_d; jy_d = x_d; jz_d = x_d;
            ns = 1;
            for i = 1:length(t_vec)
                t = t_vec(i)/slow_scale;
                if t > break_T(ns+1)
                    ns = ns+1;
                end
                t_off = t - break_T(ns);
                pos_d = pos_coef(t_off)*coeff(1+(ns-1)*p_order:ns*p_order,:);
                x_d(i) = pos_d(1); y_d(i) = -pos_d(2); z_d(i) = -pos_d(3);
                
                vel_d = (1/slow_scale)*vel_coef(t_off)*coeff(1+(ns-1)*p_order:ns*p_order,:);
                vx_d(i) = vel_d(1); vy_d(i) = -vel_d(2); vz_d(i) = -vel_d(3);
                
                acc_d = (1/slow_scale^2)*acc_coef(t_off)*coeff(1+(ns-1)*p_order:ns*p_order,:);
                ax_d(i) = acc_d(1); ay_d(i) = -acc_d(2); az_d(i) = -acc_d(3);
                
                jer_d = (1/slow_scale^3)*jer_coef(t_off)*coeff(1+(ns-1)*p_order:ns*p_order,:);
                jx_d(i) = jer_d(1); jy_d(i) = -jer_d(2); jz_d(i) = -jer_d(3);
            end
            roll_max = max(abs(atan2(ay_d,g-az_d)));
            if (roll_max > roll_MAX)
                slow_scale = slow_scale*1.1;
            end
        end
        
        yaw_d = zeros(length(t_vec),1);
        yd_d = zeros(length(t_vec),1);
        
    otherwise
        disp('Unexpected traj type');
end

%% Mapping

%Get control & nominal state
uc_nom = zeros(length(t_vec),4); %th_dot, rd, pd, yd
xc_nom = zeros(length(t_vec),10); %x,y,z, vx, vy,vz, th, r,p, y

xc_nom(:,1:6) = [x_d,y_d,z_d,vx_d,vy_d,vz_d];
xc_nom(:,10) = yaw_d;

thrust_nom = zeros(length(t_vec),1);

for i = 1:length(t_vec)
    th_vect = [ax_d(i);ay_d(i);az_d(i)-g];
    zb_des = -th_vect/norm(th_vect);
    
    thrust = norm(th_vect);
    thrust_d = -(zb_des'*[jx_d(i);jy_d(i);jz_d(i)]);
    xc_nom(i,7) = thrust;
    thrust_nom(i) = thrust;
    
    uc_nom(i,1) = thrust_d ; %thrust_dot

    %123
    roll = atan2(ay_d(i),g-az_d(i));
    pitch = asin(-ax_d(i)/thrust);
    R = rot_matrix(roll,pitch,yaw_d(i));
    
    %321
%     yc_des = [-sin(yaw_d(i)); cos(yaw_d(i)); 0];
%     xb_des = cross(yc_des,zb_des);
%     xb_des = xb_des/norm(xb_des);
%     yb_des = cross(zb_des,xb_des);
%     R = [xb_des,yb_des,zb_des];
%     roll = atan2(R(3,2),R(3,3));
%     pitch = atan2(-R(3,1),R(3,2)*sin(roll)+R(3,3)*cos(roll));
    
    xc_nom(i,8) = roll;
    xc_nom(i,9) = pitch;
    
    om_1 = [jx_d(i);jy_d(i);jz_d(i)]'*(R*[0;1;0])/thrust;
    om_2 = -[jx_d(i);jy_d(i);jz_d(i)]'*(R*[1;0;0])/thrust;
    
    %123
    roll_d = om_1*(cos(yaw_d(i))/cos(pitch)) - om_2*(sin(yaw_d(i))/cos(pitch));
    om_3 = yd_d(i) + roll_d*sin(pitch);
    om = [om_1;om_2;om_3];
    eul_rate = R_eul([roll;pitch;yaw_d(i)])*om;
    
    %321
    %roll_d = om_1 + yd_d(i)*sin(pitch);
    %pitch_d = om_2 - yd_d(i)*sin(roll)*cos(pitch);
    %eul_rate = [roll_d;pitch_d;yd_d(i)];
    
    uc_nom(i,2:4) = eul_rate';
    
end

end