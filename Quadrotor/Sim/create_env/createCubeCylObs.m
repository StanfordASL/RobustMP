function obs = createCubeCylObs(r, h, pos, numDisc)
obs = [];
for i = 1:numDisc-1
    ang = i/numDisc*pi/2;
    l = cos(ang)*r;
    w = sin(ang)*r;
    obs = [obs; createBoxObs([-l/2 -w/2 0], [l w h])];
end
obs = obs + repmat(pos,length(obs),1);
end