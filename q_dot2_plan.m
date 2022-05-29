function [q_dot2] = q_dot2_plan(q, q_dot, a, T, method, params)

% General Description:
%   This function returns a matrix with each column representing a 
%   vector of the derivatives of the joint positions with respect to time
%   at a particular time t_i

% Parameters:
% q - matrix of joint position vectors (found previously using q_plan)
% q_dot - matrix of joint velocity vectors (found previously using q_dot_plan)
% a - matrix of velocity vectors of the tool vector over T seconds
% method - user specifies whether to calculate q_dot numerically, or using
%          analytical definition of x_dot2 -JL_dot*q = JL*q_dot2
% T-time taken to move from x0 to xf
% params=[H,l2,l3]
[~, colnum] = size(q);
q_dot2 = zeros(3, colnum);
JL = zeros(3, 3, colnum);
JL_dot = zeros(3, 3, colnum);

switch method

    case 'numerical'
        
        for i = 1:3
            q_dot2(i,:) = gradient(q_dot(i,:), T/colnum);

        end

    case 'analytical'
        
        for i=1:colnum
            [~,JL(:,:,i)] = jacobian_mat_simplify(q(:,i), params);
            [~, JL_dot(:,:,i)] = jacobian_mat_dot(q(:,i),q_dot(:,i),params);
            q_dot2(:,i) = inv(JL(:,:,i))*(a(:,i)-JL_dot(:,:,i)*q_dot(:,i));
        end

end