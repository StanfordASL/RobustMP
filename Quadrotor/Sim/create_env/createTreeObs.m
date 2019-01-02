function obs = createTreeObs(pos, leavesR, treeH, numTop, numDisc)
% good numbers, numTop  = 6, numDisc = 4
trunkH = treeH/4;
trunkR = leavesR/2;

% trunk
obs = createCubeCylObs(trunkR, trunkH, [0 0 0], numDisc);

% tree top
for i = 1:numTop-1
    r = leavesR*sin(i/numTop*pi);
    obs = [obs; createCubeCylObs(2*r, (treeH-trunkH)/(numTop-1),...
        [0 0 trunkH + ((i-1)/(numTop-1))*(treeH-trunkH)], numDisc)];
end
    
obs = obs + repmat(pos,length(obs),1);
end