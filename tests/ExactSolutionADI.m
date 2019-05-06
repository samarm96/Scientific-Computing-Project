clear
v = 1; %Indexing

% Boundaries of Problem
bx = 2*pi;
by = 2*pi;
ax = 0;
ay = 0;

gridsize = 30; %Fineness of Grid
stepx = round(gridsize);
stepy = round(gridsize);
steptime = 1000;%round(.25*stepx^2);

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
               unh(1,j,t) = sin(j*dx-t*dt);
            end
        end
        for t = 1:steptime+1
            for j = 1:1:stepy+3
               un(1,j,t) = sin(j*dx-t*dt); 
            end
        end
        

%% Upper BC

        for t = 1:steptime+1
            for j = 1:stepy+3
               un(stepy+3,j,t) = sin(j*dx-2*pi-t*dt);
            end
        end
        for t = 1:steptime+1
            for j = 1:stepy+3
               unh(stepy+3,j,t) = sin(j*dx-2*pi-t*dt);
            end
        end

%% Right BC

        for t = 1:steptime+1
            for i = 1:stepx+3
               un(i,stepx+3,t) = sin(2*pi-i*dy-t*dt);
            end
        end
        for t = 1:steptime+1
            for i = 1:stepx+3
               unh(i,stepx+3,t) = sin(2*pi-i*dy-t*dt);
            end
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
tic        
        for t = 1:steptime
            %Reinitialize g and alph
            g = ones(1,N);
            alph = ones(1,N); 
            for i = 2:stepy+2            
                for j = 2:stepx+2 
                        %Preallocate variables needed for thomas algorithm
                        f = ones(1,N).*(lamb1*un(1)-2*lamb1*un(2)+lamb1*un(3))+(-sin(j*dx-i*dy-t*dt)-sin(j*dx-i*dy-t*dt)-cos(j*dx-i*dy-t*dt));
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
                end
            end
            %Reinitialize g and alph
            g = ones(1,N);
            alph = ones(1,N);            
            for i = 2:stepy+2  
                for j = 2:stepx+2  
                        %Preallocate variables needed for thomas algorithm
                        f = ones(1,N).*(lamb1*unh(1)-2*lamb1*unh(2)+lamb1*unh(3))+(-sin(j-i-t)-sin(j-i-t)-cos(j-i-t));
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
%                         uexact(i,j,t+1) = sin(j-i-t);
%                        E(v) = un(i,j,t+1) - uexact(i,j,t+1);
%                        v = v + 1;
                end
            end


        end
toc
t = toc

%   L1 = (1/(stepx*stepy*steptime))* sum(un(:,:,:) - uexact(:,j,t)); %Absolute Error
%   L2 = sqrt((1/(stepx*stepy*steptime))* sum(un(:,:,:) - uexact(:,:,:)).^2); %Mean Square Error
%   Linf = max(abs(sum(un(:,:,:) - uexact(:,:,:))));
%   L5 = (1/(stepx*stepy*steptime)).*abs(E);
%   G = 1:length(L5);
%   semilogx(G,log(L5))
%   xlabel('log(Iteration Number)');
%   ylabel('log(Error)');
%   title('Order of Error for ADI Method')
%         