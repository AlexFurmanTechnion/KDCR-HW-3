function [Mo] = backward_recursion_M(ac,ae,q,para,alpha,omega)
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

rc = [0, -h/2, 0;
      -l2/2, 0, 0;
      0,0, -l/2];
re = [0, -h, 0;
      -l2, 0, 0;
      0,0, -l];
I{4} = zeros(3,3);
I{3} = m_d4*l^2/12*eye(3);
I{3}(3,3)=0;
I{2} = m2*l2^2/12*eye(3);
I{2}(1,1)=0;
I{1} = m1*h^2/12*eye(3);
I{1}(3,3)=0;


f(:,4) = ae(:,4)*M-M*R0{3}*[0;0;-g];
Mo(:,4) =zeros(3,1);

i=3;
for j=1:3
if i >= 1
    f(:,i) = m(i)*ac(:,i) + R{i+1}*f(:,i+1)-m(i)*transpose(R0{i})*[0;0;-g];
     Mo(:,i) =R{i+1}*Mo(:,i+1)+ cross(rc(i,:),f(:,i))' + cross(re(i,:)-rc(i,:),R{i+1}*f(:,i+1))'+I{i}*alpha(:,i)+cross(omega(:,i),I{i}*omega(:,i));
   
end
i=i-1;
end


end