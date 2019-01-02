function [nn_i, Cnn_i] = get_NN(i,nodes,N,r)

%get NN for start node
nn_i = zeros(1,N);
Cnn_i = zeros(1,N);

%compute euc distances
dists = norms(nodes(i+1:N,:)-repmat(nodes(i,:),N-i,1),2,2);
nn_i(1,i+1:N) = (dists <= r)';
prune_mat(i) = 1;

%now go through list
for j = 1:N
    if prune_mat(j)
        continue;
    else
        [~,T_ij,C_ij] = steer(nodes(i,:)',nodes(j,:)',0);
        if (C_ij <= r)
            nn_i(j) = 1; Cnn_i(j) = C_ij; Tnn_i(j) = T_ij;
            [Pnn_i{1,j},~,~] = steer(nodes(i,:)',nodes(j,:)',1);
        end
    end
end


end