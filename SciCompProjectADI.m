clear

coarsemesh = 5;
r = 1.2;
flag = 1;
v =1;
m = 1;
p = 1;
l = 1;
E = 1000;
gridsize = 5;
for gridsize = [coarsemesh*r, coarsemesh, coarsemesh*r*r];
    stepx = round(gridsize);
    stepy = round(gridsize);
    steptime = 3;%round(.25*stepx^2);
    flag = 1;
    unh = ones(stepx+1, stepy+1, steptime+1); % First Step
    un = ones(stepx+1, stepy+1, steptime+1);
    N = stepx+2;
    c = ones(1,N);
    b = ones(1,N);


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
 % second step


         gabx = bx*(bx-ax)^2;
         fabx = (bx-ax)^2 * cos((pi*(bx-ax)/(bx-ax)));

        lamb = dt/(2*dx^2);
        lamb1 = dt/(2*dy^2);

        for t = 1:steptime+1
            for j = 2:stepy+2
                g = ones(1,N);
                alph = ones(1,N);                 
                for i = 2:stepx
                   
                    alph(1) = 1 + 2 * lamb;
                    g(1) = lamb1*un(1)-2*lamb1*un(2)+lamb1*un(3);
                    c(1) = -2*lamb;
                    for v = 2:N
                        alph(v) = alph(v) - (b(v)/(alph(v-1))*c(v-1));
                        g(v) = g(v) - (b(v)/(alph(v-1))*g(v-1));
                    end
                    g(N) = gabx + (j-ay)/(by-ay) * (fabx-gabx);
                    unh(N,j-1,t) = g(N) / alph(N);
                    for k = 1:N-1
                        unh(N-k,j-1,t) = g(N-k) - c(N-k) * unh(N-k+1,j-1,t)/alph(N-k);
                    end
                end
            end
        end
        for t = 1:steptime+1
            for j = 2:stepy+2
                g = ones(1,N);
                alph = ones(1,N);                 
                
                for i = 2:stepx

                    alph(1) = 1 + 2 * lamb1;
                    g(1) = lamb*unh(i-1,j-1,t)-2*lamb*unh(i,j-1,t)+lamb*unh(i+1,j-1,t);
                    c(1) = -2*lamb1;
                    for v = 2:N
                        alph(v) = alph(v) - (b(v)/(alph(v-1))*c(v-1));
                        g(v) = g(v) - (b(v)/(alph(v-1))*g(v-1));
                    end
                    g(N) = gabx + (j-ay)/(by-ay) * (fabx-gabx);
                    un(N,j-1,t) = g(N) / alph(N);
                    for k = 1:N-1
                        un(N-k,j-1,t) = g(N-k) - c(N-k) * un(N-k+1,j-1,t)/alph(N-k);
                    end
                end
            end
        end
         a1(p) = un(3,3,steptime);

         if (p > 4) == 1
            E = (a1(p-1)-a1(p-2))/a1((p-2));
         end
        if abs(E*100) <= .01
            convergence(l) = steptime;
            answers(l) = a1(p);
            l = l + 1;
            flag = 0;
        end
        p = p + 1;
        v = v+1;
        m = m + 1;
        steptime = steptime + 1;
    end
end
