function [v] = v_plan(prof,T,n,x0,xf)
%General Description:
    %This function returns a matrix, with each column representing a the
    %velocity vector of the tool in world frame at a particular time t_i

%Parameters:
    %prof: which trajectory profile the user wishes to implement
    %[NOTE: prof accepts 'constant', 'trapezoidal', and 'polynomial']
    %T: time taken to move from x0 to xf
    %n: number of velocity vectors recorded in T seconds
    %x0: initial position of tool (3x1 vector)
    %xf: final position of tool (3x1 vector)

v = zeros(3,n);
t = linspace(0, T, n);

switch prof

    case 'constant'
        for i=1:n
            v(:,i) = (xf-x0)./T;
        end

    case 'trapezoidal'
        a = (xf-x0).*(36/(5*T^2));
        
        for i=1:floor(n/6)
            v(:,i) = t(i).*a;
        end
        for i=floor(n/6)+1:floor(5*n/6)
            v(:,i) = (T/6).*a;
        end
        for i=floor(5*n/6)+1:n
            v(:,i) = (T-t(i)).*a;
        end

    case 'polynomial'
        c5 = (6/T^5).*(xf-x0);
        c4 = (-15/T^4).*(xf-x0);
        c3 = (10/T^3).*(xf-x0);

        for i=1:n
            v(:,i) = 5.*c5.*t(i)^4+4.*c4.*t(i)^3+3.*c3.*t(i)^2;
        end

end
