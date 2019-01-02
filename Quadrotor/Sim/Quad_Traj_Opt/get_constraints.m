function A = get_constraints(N_seg,T_seg,p_order,pos_coef,vel_coef,acc_coef,jer_coef,snap_coef)
%known: p_j(0,T), j = 1,...,N_seg
%       v_1(0) = v_init
%       a_1(0) = 0
%       j_1(0) = 0
%       v_Ns(T) = 0
%       j_Ns(T) = 0
%       continuity (N_seg-1) in v,a,j,s

%pos waypoints
A_wp = repmat(T_seg(1)*0,2*N_seg,p_order*N_seg);
for ns = 1:N_seg
    A_wp(1+(ns-1)*2,1+(ns-1)*p_order:ns*p_order) = pos_coef(0);
    A_wp(2*ns,1+(ns-1)*p_order:ns*p_order) = pos_coef(T_seg(ns));
end

%init velocity,acc & jerk
A_vaj_init = [vel_coef(0),repmat(T_seg(1)*0,1,p_order*(N_seg-1));
              acc_coef(0),repmat(T_seg(1)*0,1,p_order*(N_seg-1));
              jer_coef(0),repmat(T_seg(1)*0,1,p_order*(N_seg-1))];
d_vj_init = repmat(T_seg(1)*0,3,1);

%final velocity & jerk
A_vj_end = [repmat(T_seg(1)*0,1,p_order*(N_seg-1)),vel_coef(T_seg(ns));
            repmat(T_seg(1)*0,1,p_order*(N_seg-1)),jer_coef(T_seg(ns))];
d_vj_end = repmat(T_seg(1)*0,2,1);

%vel,acc,jerk,snap continuity
A_vc = repmat(T_seg(1)*0,N_seg-1,p_order*N_seg);
A_ac = A_vc;
A_jc = A_vc;
A_sc = A_vc;
d_c = repmat(T_seg(1)*0,N_seg-1,1);
for i = 1:N_seg-1
    ns = i+1;
    A_vc(i,1+(ns-1)*p_order:ns*p_order) = vel_coef(0);
    A_vc(i,1+(ns-2)*p_order:(ns-1)*p_order) = -vel_coef(T_seg(ns-1));
    A_ac(i,1+(ns-1)*p_order:ns*p_order) = acc_coef(0);
    A_ac(i,1+(ns-2)*p_order:(ns-1)*p_order) = -acc_coef(T_seg(ns-1));
    A_jc(i,1+(ns-1)*p_order:ns*p_order) = jer_coef(0);
    A_jc(i,1+(ns-2)*p_order:(ns-1)*p_order) = -jer_coef(T_seg(ns-1));
    A_sc(i,1+(ns-1)*p_order:ns*p_order) = snap_coef(0);
    A_sc(i,1+(ns-2)*p_order:(ns-1)*p_order) = -snap_coef(T_seg(ns-1));
end

A_known = [A_wp;
    A_vaj_init;
    A_vj_end;
    A_vc;
    A_ac;
    A_jc;
    A_sc];

%unknown: 
%         v_j(T), j = 1,...,N_seg-1
%         a_j(T), j = 1,...,N_seg
%         j_j(T), j = 1,....,N_seg-1
%         s_1(0)
%         s_j(T), j = 1,...,N_seg

A_v = repmat(T_seg(1)*0,N_seg-1,p_order*N_seg);
for ns = 1:N_seg-1
    A_v(ns,1+(ns-1)*p_order:ns*p_order) = vel_coef(T_seg(ns));
end
A_a = repmat(T_seg(1)*0,N_seg,p_order*N_seg);
for ns = 1:N_seg
    A_a(ns,1+(ns-1)*p_order:ns*p_order) = acc_coef(T_seg(ns));
end
A_j = repmat(T_seg(1)*0,N_seg-1,p_order*N_seg);
for ns = 1:N_seg-1
    A_j(ns,1+(ns-1)*p_order:ns*p_order) = jer_coef(T_seg(ns));
end
A_s = repmat(T_seg(1)*0,N_seg,p_order*N_seg);
for ns = 1:N_seg
    A_s(ns,1+(ns-1)*p_order:ns*p_order) = snap_coef(T_seg(ns));
end
A_s = [repmat(T_seg(1)*0,1,p_order*N_seg);A_s];
A_s(1,1:p_order) = snap_coef(0);

A_unknown = [A_v;
            A_a;
            A_j;
            A_s];

A = [A_known;
    A_unknown];

end