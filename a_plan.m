function [a] = a_plan(prof,T,n,x0,xf)
%General Description:
    %This function returns a matrix, with each column representing a the
    %acceleration vector of the tool in world frame at a particular time t_i

%Parameters:
    %prof: which trajectory profile the user wishes to implement
    %[NOTE: prof accepts 'constant', 'trapezoidal', and 'polynomial']
    %T: time taken to move from x0 to xf
    %n: number of acceleration vectors recorded in T seconds
    %x0: initial position of tool (3x1 vector)
    %xf: final position of tool (3x1 vector)

a = zeros(3,n);
t = linspace(0, T, n);

switch prof

    case 'constant'
        for i=1:n
            a(:,i) = 0;
        end

    case 'trapezoidal'
        a1 = (xf-x0).*(36/(5*T^2));
        
        for i=1:floor(n/6)
            a(:,i) = a1;
        end
        for i=floor(n/6)+1:floor(5*n/6)
            a(:,i) = 0;
        end
        for i=floor(5*n/6)+1:n
            a(:,i) = -a1;
        end

    case 'polynomial'
        c5 = (6/T^5).*(xf-x0);
        c4 = (-15/T^4).*(xf-x0);
        c3 = (10/T^3).*(xf-x0);

        for i=1:n
            a(:,i) = 20.*c5.*t(i)^3+12.*c4.*t(i)^2+6.*c3.*t(i);
        end

end