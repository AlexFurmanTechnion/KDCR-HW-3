function [X_dot]=state_eq_q6(t,X)
% dt=(2-0)/200;
index=floor(t/0.001)+1
n=load('tau_workspace.mat');
global counter
counter
tau=[n.tau(1,index) n.tau(2,index) n.tau(3,index)];
counter = counter+1;
% t1=X(1);
% t2=X(2);
% d3=X(3);
% dt1=X(4);
% dt2=X(5);
% dd3=X(6);
% Xdot=zeros(6,1);
% Xdot(1) =X(4);
% Xdot(2)= X(5);
% Xdot(3)=X(6);

rho= 7800;
d=0.015;
l=0.2;
l2=0.2; 
h=0.1; 
g=9.81;
M=0.6;
m1=pi*d^2/4*h*rho;
m2=pi*d^2/4*l2*rho;
m_d4=pi*d^2/4*l*rho; 
param=[m1 m2 m_d4 M g h l2 l ];
% T=2;
% n=200;
% tau=tau_plan(T,n);

% q=X(1:3);
% qdt=X(4:6);

d4 = X(3);
theta2=X(2);
theta1=X(1);
theta1_dt = X(4);
theta2_dt = X(5);
d4_dt = X(6);
q_dot = [theta1_dt theta2_dt d4_dt];

m = [m1 m2 m_d4 M];
c1=cos(theta1);
c2=cos(theta2);
s1=sin(theta1);
s2=sin(theta2);
a = [0 l2 0 ];
alpha = [sym(pi)/2 -sym(pi)/2 0];
d = [h 0 d4];
theta = [theta1 theta2  0 ];
dhparams = [a' alpha' d' theta'];
for i=1:3
    A(:,:,i)=[
        cos(theta(i)), -sin(theta(i))*cos(alpha(i)), sin(theta(i))*sin(alpha(i)), a(i)*cos(theta(i));
        sin(theta(i)), cos(theta(i))*cos(alpha(i)), -cos(theta(i))*sin(alpha(i)), a(i)*sin(theta(i));
        0, sin(alpha(i)), cos(alpha(i)), d(i);
        0 0 0 1
        ];
end

A01=A(:,:,1);
A02=A01*A(:,:,2);
A03=A02*A(:,:,3);

R01=A01(1:3,1:3);
R02=A02(1:3,1:3);
R03=A03(1:3,1:3);

%rotation matrix (world frame)
A0t=[c1*c2 -s1 -c1*s2   l2*c1*c2-d4*c1*s2;
     c2*s1  c1 -s1*s2   l2*c2*s1-d4*s1*s2;
     s2      0     c2   h+d4*c2+l2*s2;
     0       0      0         1];

R0t=[c1*c2 -s1 -c1*s2 ;
     c2*s1  c1 -s1*s2 ;
     s2      0     c2  ];
JL_world=[d4*s1*s2-l2*c2*s1  -d4*c1*c2-l2*c1*s2 -c1*s2; ...
          l2*c1*c2-d4*c1*s2  -d4*c2*s1-l2*s1*s2 -s1*s2;
               0             l2*c2-d4*s2          c2];
JA_world=[0 s1 0;
          0 -c1 0;
          1  0  0];

JL_tool=R0t.'*JL_world;
JA_tool=R0t.'*JA_world;
% simplify(JL_tool);
% simplify(JA_tool);

         
% jacobian of tool frame, end mass
JL_t{4} = JL_tool;
JA_t{4} = JA_tool;
% jacobian of tool frame, d4 => d4/2
JL_t{3} =   [ 0    ,  -d4/2,    0; 
    l2*c2 - d4/2*s2,    0,    0; 
            0    ,   l2,     1];
JA_t{3} = [  s2, 0, 0;
             0, -1, 0;
              c2, 0, 0];
% jacobian of tool frame, l2 => l2/2, d4=0
JL_t{2} =R02.'* [ 0    ,  0,    0; 
         l2/2*c2,    0,    0; 
            0    ,   l2/2,     0];
JA_t{2} =R02.'*[  s2, 0, 0;
             0, -1, 0;
              c2, 0, 0]  ; 
% jacobian of tool frame, h => h/2, d4=l2=0
JL_t{1} =R01.'*[ 0    ,  0,    0; 
             0,    0,    0; 
            0    ,   0,     0];
JA_t{1} =R01.'*[  0, 0, 0;
             0, 0, 0;
              0, 0, 0];



% jacobian of world frame, end mass
JL_w{4} = JL_world;
JA_w{4} = JA_world;
% jacobian of tool frame, d4 => d4/2
JL_w{3} =   [d4/2*s1*s2-l2*c2*s1  -d4/2*c1*c2-l2*c1*s2 -c1*s2; ...
          l2*c1*c2-d4/2*c1*s2  -d4/2*c2*s1-l2*s1*s2 -s1*s2;
               0             l2*c2-d4/2*s2          c2];
JA_w{3} = [0 s1 0;
          0 -c1 0;
          1  0  0];
% jacobian of tool frame, l2 => l2/2, d4=0
JL_w{2} =   [-l2/2*c2*s1  -l2/2*c1*s2 0; ...
             l2/2*c1*c2  -l2/2*s1*s2  0;
               0             l2/2*c2  0];
JA_w{2} = [0 s1 0;
          0 -c1 0;
          1  0  0]; 
% jacobian of tool frame, h => h/2, d4=l2=0
JL_w{1} =  [0  0  0; ...
            0  0  0;
            0  0  0];
JA_w{1} = [0 0 0;
          0  0 0;
          1  0  0];

% moment of interia
I{4} = zeros(3,3);
I{3} = m_d4*l^2/12*eye(3);
I{3}(3,3)=0;
I{2} = m2*l2^2/12*eye(3);
I{2}(1,1)=0;
I{1} = m1*h^2/12*eye(3);
I{1}(3,3)=0;

H = zeros(3,3);
for i = 1:4
    H = H + m(i)*JL_t{i}.'*JL_t{i} + JA_t{i}.'*I{i}*JA_t{i};
end
% simplify(H)
G = zeros(3,1);
for i=1:4
    G = G - m(i)*JL_w{i}.'*[0;-g;0];
end
% simplify(G)

C = zeros(3,3);

% Christoff = zeros(3,3);
Chr(:,:,1) = zeros(3,3);
Chr(:,:,2) = [M*d4^2*sin(2*theta2)-M*l2^2*sin(2*theta2)+(d4^2*m_d4*sin(2*theta2))/3-(l2^2*m2*sin(2*theta2))/3-l2^2*m_d4*sin(2*theta2)-2*M*d4*l2*cos(2*theta2)-d4*l2*m_d4*cos(2*theta2)    0     0;
    0 0 0;
    0 0 0];
Chr(:,:,3) = [M*d4+(d4*m_d4)/3-M*d4*cos(2*theta2)-M*l2*sin(2*theta2)-(d4*m_d4*cos(2*theta2))/3-(l2*m_d4*sin(2*theta2))/2  0  0;
    0    2*M*d4+2*d4*m_d4/2  0;
    0   0   0];

for i=1:3
    for j=1:3
        for k = 1:3
            C(i,j) = C(i,j)+ 1/2*(Chr(i,j,k)+Chr(i,k,j)-Chr(j,k,i))*q_dot(k);
        end
    end
end
% simplify(C)


PIH=zeros(1,2);
qddt=inv(H)*(tau'-C*q_dot'-G);
% q_dot
% qddt
% size(qddt)
% tauu=tau
% q_dot
% qddt
% gee=G
% mult=C*q_dot'
% ivH=inv(H)
% tau-C*q_dot'-G
X_dot=double([q_dot,qddt']')

end