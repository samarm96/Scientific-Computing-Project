clear

% Boundaries of Problem
bx = 2*pi;
by = 2*pi;
ax = 0;
ay = 0;

gridsize = 25; %Fineness of Grid
stepx = gridsize;
stepy = gridsize;
steptime = 1000; %Number of timesteps

unh = zeros(stepx+3, stepy+3, steptime+1); % First Step
un = zeros(stepx+3, stepy+3, steptime+1); %Second Step
N = stepx+2;

%Find step size
dx = (bx-ax)/(stepx+2);          
dy = (by-ay)/(stepy+2);          
dt = dx^2/4;                
%% Bottom BC

        for t = 1:steptime+1
            for j = 1:1:stepy+3
               unh(1,j,t) = j.*(j-ax).^2; 
            end
        end
        for t = 1:steptime+1
            for j = 1:1:stepy+3
               un(1,j,t) = j.*(j-ax).^2; 
            end
        end
        

%% Upper BC

        for t = 1:steptime+1
            for j = 1:stepy+3
               un(stepy+3,j,t) = (i-ax)^2 * cos((pi*(j-ax)/(bx-ax)));
            end
        end
        for t = 1:steptime+1
            for j = 1:stepy+3
               unh(stepy+3,j,t) = (i-ax)^2 * cos((pi*(j-ax)/(bx-ax)));
            end
        end

%% Right BC
        gabx = bx*(bx-ax)^2;
        fabx = (bx-ax)^2 * cos((pi*(bx-ax)/(bx-ax)));

        for t = 1:steptime+1
            for i = 1:stepx+3
               un(i,stepx+3,t) = gabx + (i-ay)/(by-ay) * (fabx - gabx);
            end
        end
        for t = 1:steptime+1
            for i = 1:stepx+3
               unh(i,stepx+3,t) = gabx + (i-ay)/(by-ay) * (fabx - gabx);
            end
        end
        
if exist( 'checkpoint.mat','file' ) % If a checkpoint file exists, load it
    fprintf('Checkpoint file found - Loading\n');
    load('checkpoint.mat')
end        
%% ADI Implicit Calcs

        %Preallocate variables needed for thomas algorithm that don't need
        %to be reinitialized
        lamb = dt/(2*dx^2);
        lamb1 = dt/(2*dy^2);
        c = ones(1,N).*lamb;
        c(1) = -2*lamb;
        c(N) = 0;        
        b = ones(1,N).*lamb;
        b(1) = 0;
        
        for t = 1:steptime
            %Reinitialize g and alph
            g = ones(1,N);
            alph = ones(1,N); 
            for i = 2:stepy+1            
                for j = 2:stepx+2 
                        %Preallocate variables needed for thomas algorithm
                        f = ones(1,N).*(lamb1*un(1)-2*lamb1*un(2)+lamb1*un(3));
                        g(1) = f(1);                   
                        a = ones(1,N) .* 1+2*lamb;
                        a(1) = 1 + 2 * lamb;
                        alph(1) = a(1);
                        a(N) = -2*lamb;
                        
                        %Thomas Algorithm
                        for v = 2:N
                            alph(v) = a(v) - (b(v)/(alph(v-1))*c(v-1));
                            g(v) = f(v) - (b(v)/(alph(v-1))*g(v-1));
                        end
                        unh(i-1,N,t) = g(N) / alph(N);
                        for k = 1:N-1
                            unh(i-1,N-k,t) = (g(N-k) - c(N-k) * unh(i,N-k+1,t))/alph(N-k);
                        end
                    for r = 2:stepx+1             % Neumann BC
                        unh(stepy+2,r,t) = unh(stepy+1,r,t);
                    end                          
                end
            end
            %Reinitialize g and alph 
            g = ones(1,N);
            alph = ones(1,N);            
            for i = 2:stepy+1  
                for j = 2:stepx+2  
                        %Preallocate variables needed for thomas algorithm                
                        f = ones(1,N).*(lamb1*unh(1)-2*lamb1*unh(2)+lamb1*unh(3));
                        g(1) = f(1);
                        a = ones(1,N) .* 1+2*lamb;
                        a(1) = 1 + 2 * lamb;
                        alph(1) = a(1);
                        a(N) = -2*lamb;
                        
                        %Thomas Algorithm
                        for v = 2:N
                            alph(v) = a(v) - (b(v)/(alph(v-1))*c(v-1));
                            g(v) = f(v) - (b(v)/(alph(v-1))*g(v-1));
                        end
                        un(i-1,N,t) = g(N) / alph(N);
                        for k = 1:N-1
                            un(i-1,N-k,t) = (g(N-k) - c(N-k) * un(i,N-k+1,t))/alph(N-k);
                        end
                end
                for r = 2:stepx+1             % Neumann BC
                    un(stepy+2,r,t) = un(stepy+1,r,t);
                end                  
            end
            
    if mod(t,t/4)==0
        %save checkpoint   
        fprintf('Saving checkpoint\n');
        save('checkpoint.mat');
    end
        end

