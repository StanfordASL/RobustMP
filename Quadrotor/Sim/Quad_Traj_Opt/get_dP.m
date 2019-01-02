function d_P = get_dP(d_F,s_F,T_seg,Q,A_inv)

Q_full = Q(T_seg);
R = A_inv'*Q_full*A_inv;
R_PP = R(s_F+1:end,s_F+1:end);
R_FP = R(1:s_F,s_F+1:end);
d_P = -R_PP\(R_FP'*d_F);

end