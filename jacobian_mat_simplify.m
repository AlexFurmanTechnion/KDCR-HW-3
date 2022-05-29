
function [J, JL] = jacobian_mat_simplify(qi,params) 
% jacobian_mat computes jacobian matrix of the robot given configuration q
%q - 3x1 vector of current joint values
% params=[H,l2,l3]
l2 = params(2);
l3 = params(3);

c1 = cos(qi(1));
s1 = sin(qi(1));
c2 = cos(qi(2));
s2 = sin(qi(2));
d4=qi(3);
sig1=c1*c2*l2-c1*s2*(d4+l3);
sig2=s2*l2+c2*(d4+l3);
sig3=s1*s2*(d4+l3);
sig4=c2*s1*l2;
J =     [sig3-sig4 , -(c1)*sig2 , -c1*s2 ;
        sig1 , -(s1)*sig2 , -s1*s2 ;
        0 , c1*sig1+s1*(sig4-sig3) , c2 ;
        0 , s1 , 0  ;
        0 , -c1 , 0;
        1 , 0 , 0 ];

JL = J(1:3,:);
end 