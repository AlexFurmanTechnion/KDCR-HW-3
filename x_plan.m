function [x] = x_plan(prof,T,n,x0,xf)
%General Description:
    %This function returns a matrix, with each column representing the
    %position vector of the tool in world frame at a particular time t_i

%Parameters:
    %prof: which trajectory profile the user wishes to implement
    %[NOTE: prof accepts 'constant', 'trapezoidal', and 'polynomial']
    %T: time taken to move from x0 to xf
    %n: number of position vectors recorded in T seconds
    %x0: initial position of tool (3x1 vector)
    %xf: final position of tool (3x1 vector)

x = zeros(3,n);
t = linspace(0, T, n);

switch prof

    case 'constant'
        for i=1:n
            x(:,i) = x0 + ((xf-x0)./T).*t(i);
        end

    case 'trapezoidal'
        a = (xf-x0).*(36/(5*T^2));
        
        for i=1:floor(n/6)
            x(:,i) = ((1/2)*(t(i)^2)).*a + x0;
        end
        for i=floor(n/6)+1:floor(5*n/6)
            x(:,i) = ((T/6)*t(i)-(T^2)/72).*a + x0;
        end
        for i=floor(5*n/6)+1:n
            x(:,i) = (-(1/2)*t(i)^2+T*t(i)-(13/36)*T^2).*a + x0;
        end

    case 'polynomial'
        c5 = (6/T^5).*(xf-x0);
        c4 = (-15/T^4).*(xf-x0);
        c3 = (10/T^3).*(xf-x0);

        for i=1:n
            x(:,i) = c5.*t(i)^5+c4.*t(i)^4+c3.*t(i)^3+x0;
        end
end

end
