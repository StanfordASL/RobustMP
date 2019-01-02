%print matrices for copying over to ROS C++ (Metric.h)

clear; close all; clc;

load 'Quad_Metric_Coeffs.mat';

n_mon = size(W_sol,3);

%% Coeffs
for k = 1:n_mon
    fprintf('(Eigen::Matrix<double,9,9>() <<');
    for i = 1:9
        fprintf('%.8f,',W_sol(i,:,k)); 
        if (i<9)
            fprintf('\n');
        else
            fprintf('\b');
        end
    end
    fprintf(').finished(),\n');
end

%% Powers

p = full(p); p = 1.0*p;
for i = 1:n_mon
    fprintf('%.2f,',p(i,:)); fprintf('\n');
end

