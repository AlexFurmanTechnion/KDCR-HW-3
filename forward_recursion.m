function [omega,alpha,ae,ac] = forward_recursion(q, q_dot, q_ddot,para)
%para=[m1 m2 m_d4 M g h l2 l ];
m1=para(1);
m2=para(2);
m_d4=para(3);
M=para(4);
g=para(5);
h=para(6);
l2=para(7);
l=para(8);

d4 = q(3);
theta2=q(2);
theta1=q(1);
m = [m1 m2 m_d4 M];

A01 = [cos(theta1) 0 sin(theta1) 0;
    sin(theta1)    0 -cos(theta1) 0;
    0 1 0 h;
    0 0 0 1];
A12 = [cos(theta2) 0 -sin(theta2) l2*cos(theta2); 
       sin(theta2) 0 cos(theta2)  l2*sin(theta2);
        0 -1 0 0;
        0 0 0 1];
A23 = [1 0 0 0;
       0 1 0 0;
       0 0 1 l;
       0 0 0 1];

R{1} = A01(1:3,1:3);
R{2} = A12(1:3,1:3);
R{3} = A23(1:3,1:3);
u = [0,0,0;
    1,-1,0;
    0,0,1];

theta_dot = [q_dot(1),q_dot(2),0];
d_dot = [0,0,q_dot(3)];

theta_ddot = [q_ddot(1),q_ddot(2),0];
d_ddot = [0,0,q_ddot(3)];

rc = [0, -h/2, 0;
      -l2/2, 0, 0;
      0,0, -l/2];
re = [0, -h, 0;
      -l2, 0, 0;
      0,0, -l];

    omega  = zeros(3,4); 
    alpha  = zeros(3,4);
    ac     = zeros (3,4);
    ae     = zeros(3,4);
    for i=1:1:3
    omega(:,i+1) = transpose(R{i})*omega(:,i) + u(:,i)*theta_dot(i);
    alpha(:,i+1) = transpose(R{i})*alpha(:,i) + u(:,i)*theta_ddot(i) + cross(omega(:,i),u(:,i)*theta_dot(i));
    ac(:,i+1) = transpose(R{i})*ae(:,i) + cross(alpha(:,i), rc(i,:))' + cross(omega(:,i), cross(omega(:,i),rc(i,:)))' +...
    u(:,i)*d_ddot(i) + 2*cross(omega(:,i),u(:,i))*d_dot(i);
    ae(:,i+1) = transpose(R{i})*ae(:,i) + cross(alpha(:,i), re(i,:))' + cross(omega(:,i), cross(omega(:,i),re(i,:)))' +...
    u(:,i)*d_ddot(i) + 2*cross(omega(:,i),u(:,i))*d_dot(i);
    end

end