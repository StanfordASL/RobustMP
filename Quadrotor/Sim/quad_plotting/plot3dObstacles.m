function plot3dObstacles(obstacles, varargin)
color = [1 0.6 0.6];
alpha = 0.7;
enhance = 0;

if nargin > 3
    color = varargin{1};
    alpha = varargin{2};
    enhance = varargin{3};
elseif nargin > 2
    color = varargin{1};
    alpha = varargin{2};
elseif nargin == 2
    color = varargin{1};
end
    
for i=1:2:length(obstacles(:,1))
    a = obstacles(i,1);
    b = obstacles(i,2);
    c = obstacles(i,3);
    d = obstacles(i+1,1);
    e = obstacles(i+1,2);
    f = obstacles(i+1,3);
    
    plotcube([a-d b-e c-f],[d e f],alpha,color,enhance);
end
end