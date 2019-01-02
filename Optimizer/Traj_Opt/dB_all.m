function dB_all = dB_all(dB,x,u,n,m,N)

dB_all = zeros((N+1)*n);

for j = 1:N+1
   xj = x(1+(j-1)*n:j*n);
   uj = u(1+(j-1)*m:j*m);
   dB_j = dB(xj);
   dB_all(1+(j-1)*n:j*n,1+(j-1)*n:j*n) = sum(dB_j.*reshape(uj,1,1,m),3);
end

end