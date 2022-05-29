function [f]=backward_recursion(ac,ae,q,para)
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
R{4} = eye(3);

R0{1} = R{1};
R0{2} = R0{1}*R{2};
R0{3} = R0{2}*R{3};

f(:,4) = ae(:,4)*M-M*R0{3}*[0;0;-g];

i=3;
if i >= 1
    f(:,i) = m(i)*ac(:,i) + R{i+1}*f(:,i+1)-m(i)*transpose(R0{i})*[0;0;-g];
    i=i-1;
end

end