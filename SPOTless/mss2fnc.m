function fn=mss2fnc(p,x,v)
% function q=dmsubs(p,x,v)
%
%
% Data matrix substitution.
%
%
% INPUTS:
%   p  -  m-by-1 msspoly
%   x  -  k-by-1 free msspoly
%   v  -  k-by-n real double
%
% OUTPUT:
%   q  -  m-by-n double
%
% DESCRIPTION: q(:,i) is the result of substituting v(:,i) for x in p
% p must be an msspoly in x alone.

if nargin < 3, error('Three arguments required.'); end
if ~isa(p,'msspoly'), error('First argument must be an msspoly.'); end

m = size(p,1);
if size(p,2) ~= 1, error('First argument must be m-by-1.'); end

[f,xn] = isfree(x);
if ~f, error('Second argument must be free msspoly.'); end
k = size(x,1);
if size(x,2) ~= 1, error('Second argument must by k-by-1'); end

if size(v,1) ~= k,
    error(['Second/Third arguments must have the same number of rows']);
end

if ~isa(v,'double'), error('Third argument must be a double.'); end

% fn_vec = [];
str = '@(x) [';
for idx = 1:m
    
    [xd,pd,Md] = decomp(p(idx));
    
    
    N = size(v,2);
    po = size(pd,1);
    n = size(pd,2);
    
    if n == 0
        q = repmat(double(p(idx)),1,N);
        str = [str sprintf('%f*ones(1,size(x,2));',q(1))];
        continue;
    end
    
    % First, test that xd is a subset of xn.
    % Sort out indicies.
    [~,xdn] = isfree(xd);
    perm = mss_match(xn,xdn);
    perm = perm(perm ~= 0);
    
    if length(perm) ~= length(xd)
        error('p must only be a function of x');
    end
    
    % [~,I] = sort(perm);
    % pd = pd(:,I);
    
    % str = '@(coef,x) [';
    
    for k=1:size(Md,1),
        for i=1:size(Md,2),
            if i > 1,
                str = [str '+'];
            end
            %     str = [str sprintf('coef(%d,%d)',k,i)];
            str = [str sprintf('%f',full(Md(k,i)))];
            for j=1:size(pd,2),
                if pd(i,j) ~= 0
                    str = [str sprintf('.*x(%d,:).^%d',perm(j),full(pd(i,j)))];
                end
            end
        end
        str = [str ';'];
    end
end
    str = [str ']'];
    fn = str2func(str);
end
