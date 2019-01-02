function GradObj = Geodesic_grad_MPC(~,w,K,N,n,Ti,W,dW,~,Phi_dot,n_W)

global GEO_X_MPC;
global GEO_MXDOT_MPC;

GradObj = zeros(n*(N+1),1);

dW_poly = dW(GEO_X_MPC);

for k = 1:K+1 %0 ---> K
    M_xdot = GEO_MXDOT_MPC(:,k);
    
    GradObj = GradObj + w(k)*Phi_dot(:,:,k)'*M_xdot;
    
%     W_dx = dW(GEO_X_MPC(:,k));
    
    for j = 1:length(n_W)
        GradObj = GradObj+...
           -(1/2)*w(k)*(M_xdot'*W.W_eval(dW_poly{j}(:,k))*M_xdot)*Ti(:,k,j);
    end
    
end


return