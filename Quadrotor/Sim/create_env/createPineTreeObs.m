function obs = createPineTreeObs(pos, leavesR, treeH, numTop, numDisc)
% good values, numTop = 6, numDisc = 4
trunkH = treeH/4;
trunkR = leavesR/4;

% trunk
obs = createCubeCylObs(2*trunkR, trunkH, [0 0 0], numDisc);

% tree top
for i = 0:numTop-1
    obs = [obs; createCubeCylObs(2*(1-i/numTop)*leavesR, (treeH-trunkH)/numTop,...
        [0 0 trunkH + (i/numTop)*(treeH-trunkH)], numDisc)];
end
    
obs = obs + repmat(pos,length(obs),1);
end