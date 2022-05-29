function [q] = q_plan(x, elbows, params)

% General Description: 
%     This function returns a matrix with each column representing the 
%     joint position vector at a particular time t_i

% Parameters:
%     x: matrix of preplanned position vectors, forming a trajectory
%     elbows: 1x2 matrix which decides

[~,colnum]=size(x);
q = zeros(3,colnum);

for i=1:colnum
    q(:,i) = inverse_kin(x(:,i), elbows, params);
end
