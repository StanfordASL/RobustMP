function ceq = brach_con(xt,D,N,n)

Y = xt(2:2:2*(N+1));
theta = xt(2*(N+1)+1:end-1);

% c = [];
ceq = (2/xt(end))*D*xt(1:2*(N+1)) - ...
    sqrt(2)*(kron(sqrt(Y).*cos(theta),[1;0])+...
             kron(sqrt(Y).*sin(theta),[0;1]));
                                           

end