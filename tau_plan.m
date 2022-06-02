function [tau]=tau_plan(T,n)

global l2 l3 h m1 m2 m_d4 M g rho d l
rho= 7800;
d=0.015;
l=0.2;
l2=0.2; 
h=0.1; 
g=9.81;
M=0.5;
m1=pi*d^2/4*h*rho;
m2=pi*d^2/4*l2*rho;
m_d4=pi*d^2/4*l*rho; 
para=[m1 m2 m_d4 M g h l2 l ];

PROFILE = 'constant';
T = T;
n = n;
h = 0.1;
l2 = 0.2;
l3 = 0;
params = [h,l2,l3];
x_init = [0.1, 0.2, 0.35]';
x_fin = [0.3, -0.1, 0.15]';
x = x_plan(PROFILE, T, n, x_init, x_fin);
v = v_plan(PROFILE, T, n, x_init, x_fin);
a = a_plan(PROFILE, T, n, x_init, x_fin);

q = q_plan(x, [1,1], params);

q_dot = q_dot_plan(q, v, T, 'analytical', params);

q_dot2 = q_dot2_plan(q, q_dot, a, T,'analytical', params);

[H,C,G] = dynamics_mat(q,q_dot,para);
H=round(H,4)
C=round(C,4)
G=round(G,4)
for i=1:200
tau(:,i) = H(:,:)*q_dot2(:,i) + C(:,:)*q_dot(:,i) + G;
end
end