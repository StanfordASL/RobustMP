function [X_e, X_dot_e,J_opt,exitflag,Result,Prob] = compute_geodesic(Prob,n,N,start_p,end_p,T_e,T_dot_e,warm,solver)

%Compute geodesic

%% update endpoints

beq = [start_p; end_p];

%if this line gives problems, comment out the crossover check
Prob = modify_b_L(Prob,beq,1:2*n);
Prob = modify_b_U(Prob,beq,1:2*n);

%% recall warm-start solution

if (~warm.sol)
    %construct straight line initial guess
    vars_0 = zeros(n*(N+1),1);
    for i = 1:n
        vars_0((i-1)*(N+1)+1:(i-1)*(N+1)+2) = [(start_p(i)+end_p(i))/2;
            -(start_p(i)-end_p(i))/2];
    end
    Prob = modify_x_0(Prob,vars_0);
else
    %warm start
    Prob = WarmDefSOL(solver,Prob,warm.result);
end

%% get ready for solve

if ~Prob.CHECK
    Prob = ProbCheck(Prob,solver);
end

%% Solve

if strcmp(solver,'npsol')
    Result = npsolTL(Prob);
else
    Result = snoptTL(Prob);
end

%% Parse

C_opt = (reshape(Result.x_k,N+1,n))';
X_e = C_opt*T_e;
X_dot_e = 2*C_opt*T_dot_e;

J_opt = Result.f_k;
exitflag = Result.Inform;%GOOD: {0,1,6}

end



