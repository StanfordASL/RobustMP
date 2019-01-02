%% Collision Checker
% Written by Brian Ichter
% Given a v and end state, check whether the resulting edge is free of
% collisions
%
% Given:
%   v- start state
%   w- end state
%   obstacles- obstacle set
%
% Return:
%   valid- true if the path between v and w is valid
%   freeDist- distance to obstacle path is invalid
%

function [valid] = checkCollision(v, w, obstacles)
bb_min = min(v,w);
bb_max = max(v,w);

for k = 1:2:length(obstacles(:,1))
    obs = obstacles(k:k+1,:);
    if ~BroadphaseValidQ(bb_min, bb_max,obs)
        if ~MotionValidQ(v, w, obs)
            valid = 0;
            return;
        end
    end
end
valid = 1;

    function [isValid] = BroadphaseValidQ(bb_min, bb_max, obs)
        for i = 1:size(obs,2)
            if bb_max(i) <= obs(1,i)  || obs(2,i) <= bb_min(i)
                isValid = 1;
                return
            end
        end
        isValid = 0;
    end

    function [isValid] = MotionValidQ(v, w, obs)
        v_to_w = w-v;
        corner = v < obs(1,:);
        lambdas = (corner.*obs(1,:) + ~corner.*obs(2,:) - v)./v_to_w+5*eps;
        for i = 1:size(obs,2)
            if FaceContainsProjection(v, v_to_w, lambdas(i), i, obs)
                isValid = 0;
                minLambda = min(abs(lambdas));
                return
            end
        end
        isValid = 1;
    end

    function [isValid] = FaceContainsProjection(v, v_to_w, lambda, j, obs)
        for i = 1:size(obs,2)
            if i ~= j && ~(obs(1,i) <= v(i) + v_to_w(i)*lambda && ...
                    v(i) + v_to_w(i)*lambda <= obs(2,i))
                isValid = 0;
                return
            end
        end
        isValid = 1;
    end
end