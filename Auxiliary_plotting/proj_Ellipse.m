function h = proj_Ellipse(dim, A, bound,center, N, c, varargin)

if nargin > 6
    th = varargin{1};
else
    th = 0;
end
n = size(A,1);

In = eye(n);
P = In(dim,:);
Ap = (P*(A\P'))\eye(length(dim));
h = Ellipse_plot(Ap*(1/bound),center,N,c,0.8,th);

end