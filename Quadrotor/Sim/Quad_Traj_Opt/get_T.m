function T_seg_up = get_T(G_T_fun,k_T,T_seg,d_F,d_P,n_T,step_T)

N_seg = length(T_seg);
dT = zeros(N_seg,1);

d = [d_F;d_P];

for it = 1:n_T
    
    %get grad
    for ns = 1:N_seg
        dT(ns) = trace(d'*G_T_fun{ns}(T_seg)*d) + k_T;
    end
    
    %take step
    T_seg = max(T_seg - step_T*dT,0.1);
    
end

T_seg_up = T_seg;

end