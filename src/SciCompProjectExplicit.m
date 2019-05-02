%% Sci Computing Project

clear
v = 1;
r = 2;
coarsemesh = 5;
a1 = ones(1,3);
flag = 1;
dtime = 0;
E = 1000;
m = 1;
l = 1;
convergence = zeros(1,3)

   for gridsize = [coarsemesh*r*r, coarsemesh*r, coarsemesh];
        stepx = gridsize;
        stepy = gridsize;
        steptime = round(.25*stepx^2);
        flag = 1;
       while flag == 1

        bx = 2*pi;
        by = 2*pi;
        ax = 0;
        ay = 0;

        %fa = (x-ax)^2 * cos((pi(x-ax)/(bx-ax)));
        %ga = x(x-ax)^2;
        t = 1;
        dx = (bx-ax)/(stepx+1);          
        dy = (by-ay)/(stepy+1);          
        dt = t/(steptime+1);          
        u = zeros(stepx+1, stepy+1, steptime+1);

        %% bottom boundary condition
        % = g(a) at the bottom

        %u(i,j,t)
        for t = 1:steptime+1
            for i = 1:1:stepy+1
               u(1,i,t) = i*(i-ax)^2; 
            end
        end
        %% Upper Boundary Condition
        % = f(a)
        for t = 1:steptime+1
            for i = 1:stepy+1
               u(stepy+1,i,t) = (i-ax)^2 * cos((pi*(i-ax)/(bx-ax)));
            end
        end


        %% Right Boundary Condition
        gabx = bx*(bx-ax)^2;
        fabx = (bx-ax)^2 * cos((pi*(bx-ax)/(bx-ax)));

        for t = 1:steptime+1
            for j = 1:stepx+1
               u(j,stepx+1,t) = gabx + (j-ay)/(by-ay) * (fabx - gabx);
            end
        end
        %% Explicit Calculation

        for t = 1:1:steptime
            for j = 2:1:stepx
                for i = 2:1:stepy
                    u(j,i,t+1) = (dt/dx^2)*(u(i,j+1,t) - 2*u(i,j,t) + u(i,j-1,t)) + (dt/dy^2)*(u(i+1,j,t) - 2*u(i,j,t) + u(i-1,j,t));  
                end
            end
            for i = 2:stepy
                u(j+1,1,t) = 0;
            end
        end
        a1(v) = u(3,3,steptime);

         if (v > 4) == 1
            E = (a1(v-1)-a1(v-2))/a1((v-2));
         end
        if abs(E*100) <= .01
            convergence(l) = steptime;
            l = l + 1;
            flag = 0;
        end
        
        v = v+1;
        m = m + 1;
        steptime = steptime + 1;
       end
    
    end
p = log((a1(3)-a1(2))/(a1(2)-a1(1))/log(r));
fexact = a1(1) - (a1(2)-a1(3))/(r^p-1);

