function [N_nn,C_nn] = compute_NN(nodes,r)

%compute cost-threshold nearest neighbors
%ds: step-size
%r: connection radius
%nodes: N x 3

%output:
%N_nn: nearest-neighbor connection graph
%C_nn: associated costs (euclidean distances)

N = size(nodes,1);

%% Compute matrix of euclidean norm differences

N_nn = zeros(N); 
C_nn = zeros(N);
for i = 1:N
   dists = norms(nodes(i+1:N,:)-repmat(nodes(i,:),N-i,1),2,2);
   N_nn(i,i+1:N) = (dists <= r)';
   N_nn(i+1:N,i) = N_nn(i,i+1:N)';
   
   C_nn(i,i+1:N) = N_nn(i,i+1:N).*dists';
   C_nn(i+1:N,i) = C_nn(i,i+1:N)';
end

N_nn = sparse(N_nn);
C_nn = sparse(C_nn);

end

