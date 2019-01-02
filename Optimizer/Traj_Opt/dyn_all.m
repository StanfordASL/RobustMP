function dyn_all = dyn_all(f,B,x,u,n,m,N)

global B_full;
B_full = kron(eye(N+1),zeros(n,m));

f_all = zeros((N+1)*n,1);
for j = 1:N+1
    xj = x(1+(j-1)*n:j*n);
    f_all(1+(j-1)*n:j*n) = f(xj);
    B_full(1+(j-1)*n:j*n,1+(j-1)*m:j*m) = B(xj);
end

dyn_all = f_all + B_full*u;

end