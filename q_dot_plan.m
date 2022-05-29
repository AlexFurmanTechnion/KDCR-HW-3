function [q_dot] = q_dot_plan(q, v, T, method, params)

% General Description:
%   This function returns a matrix with each column representing a 
%   vector of the derivatives of the joint positions with respect to time
%   at a particular time t_i

% Parameters:
% q - matrix of joint position vectors (found previously using q_plan)
% v - matrix of velocity vectors of the tool vector over T seconds
% method - user specifies whether to calculate q_dot numerically, or using
%          analytical definition of x_dot = JL*q_dot
% T-time taken to move from x0 to xf
% params=[H,l2,l3]

[~, colnum] = size(q);
q_dot = zeros(3, colnum);
JL = zeros(3, 3, colnum);
switch method

    case 'numerical'
        
        for i = 1:3
            q_dot(i,:) = gradient(q(i,:), T/colnum);

        end

    case 'analytical'
        
        for i=1:colnum
            [~,JL(:,:,i)] = jacobian_mat_simplify(q(:,i), params);
            q_dot(:,i) = inv((JL(:,:,i)))*v(:,i);
        end


end
