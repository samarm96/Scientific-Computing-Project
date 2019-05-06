%% Sci Computing Project

clear

% Boundaries of Problem
bx = 2*pi; %2pi
by = 2*pi; %2pi
ax = 0;
ay = 0;


gridsize = 10; %Fineness of Grid
stepx = round(gridsize);
stepy = round(gridsize);
steptime = 100; %Number of Timesteps

u = zeros(stepx+3, stepy+3, steptime+1); %Preallocate u matrix
uexact = zeros(stepx+3, stepy+3, steptime+1);  %Preallocate uexact matrix

v = 1; %Indexing purposes

%Find step size
dx = (bx-ax)/(stepx+2);          
dy = (by-ay)/(stepy+2);          
dt = dx^2/4;

%Define for surface plotting
x= 0:dx:(bx-ax); 
y= 0:dy:(by-ay);
x = repmat(x,1,1,steptime);
y = repmat(y,1,1,steptime);        
%% Bottom BC

        for t = 1:steptime+1
            for j = 2:1:stepx+2
               u(1,j,t) = sin(j*dx-t*dt);
            end
        end
        

%% Upper BC
        
        for t = 1:steptime+1
            for j = 2:stepx+3
               u(stepy+2,j,t) = sin(j*dx-2*pi-t*dt);
            end
        end

%% Left BC
        for t = 1:steptime+1
            for i = 2:stepy+3
               u(i,stepx+2,t) = sin(-i*dx-t*dt);
            end
        end
%% Right BC

        for t = 1:steptime+1
            for i = 2:stepy+3
               u(i,stepx+2,t) = sin(2*pi-i*dy-t*dt);
            end
        end
%% Explicit Calc
        for t = 1:1:steptime+1
            for i = 2:1:stepy+2
              for j = 2:1:stepx+2
                           
                    u(i,j,t+1) = (u(i,j,t) + (dt/dx^2)*(u(i,j+1,t) - 2*u(i,j,t) + u(i,j-1,t)) + (dt/dy^2)*(u(i+1,j,t) - 2*u(i,j,t) + u(i-1,j,t)))+(-sin(j*dx-i*dy-t*dt)-sin(j*dx-i*dy-t*dt)-cos(j*dx-i*dy-t*dt));   
                    uexact(i,j,t+1) = sin(j*dx-i*dy-t*dt); %Exact U (from manufactured solution)
                    E(v) = u(i,j,t+1) - uexact(i,j,t+1); % Difference between u and uexact for each step
                    v = v + 1; %Indexing

              end
            end
            %Used to plot the numerical solution
%             h = surf(x(:,:,1),y(:,:,1),u(:,:,t),'EdgeColor','none');
%             colorbar
%             shading faceted
%             axis ([0 2*pi 0 2*pi -180 70])
%             title({['time = ',num2str(t*dt)]})
%             xlabel('x')
%             ylabel(' y')
%             zlabel('Numerical Solution to the Manufactured Solution  \rightarrow')
%             drawnow; 
%             refreshdata(h)
        end
    
        %Used to find error
%   L1 = (1/(stepx*stepy*steptime))* sum(u(:,:,:) - uexact(:,j,t)); %Absolute Error
%   L2 = sqrt((1/(stepx*stepy*steptime))* sum(u(:,:,:) - uexact(:,:,:)).^2); %Mean Square Error
%   Linf = max(abs(sum(u(:,:,:) - uexact(:,:,:))));
%   L5 = sqrt((1/(stepx*stepy*steptime)).*abs(E).^2);
%   G = 1:length(L5);
%   semilogx(G,log(L5))
%   xlabel('log(Iteration Number)');
%   ylabel('log(Error)');
%   title('Order of Error for Explicit Method')


    

