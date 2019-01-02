function obs_coll = corner2Hull(obs,n_obs)
%input: ll/ur format
%output: convex hull

obs_coll = zeros(8,3,n_obs);

for k = 1:n_obs
    ll = 1+(k-1)*2;
    ur = ll+1;
    
    l = obs(ur,1)-obs(ll,1);
    w = obs(ur,2)-obs(ll,2);
    h = obs(ur,3)-obs(ll,3);
    
    obs_coll(:,:,k) = [obs(ll,:);
                       obs(ll,:)+[l,0,0];
                       obs(ll,:)+[l,w,0];
                       obs(ll,:)+[0,w,0];
                       obs(ll,:)+[0,0,h];
                       obs(ll,:)+[l,0,h];
                       obs(ll,:)+[l,w,h];
                       obs(ll,:)+[0,w,h]];
                       
end


end