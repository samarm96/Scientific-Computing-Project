%% Sci Computing Project

clear


% Boundaries of Problem
bx = 2*pi; 
by = 2*pi; 
ax = 0;
ay = 0;

gridsize = 10; %Fineness of grid
stepx = round(gridsize);
stepy = round(gridsize);
steptime = 300; %Number of timesteps

u = zeros(stepy+3, stepx+3, steptime+1); %Preallocate u matrix

%Find step size
dx = (bx-ax)/(stepx+2);          
dy = (by-ay)/(stepy+2);          
dt = dx^2/4;  

%Define for surface plotting
x= 0:dx:(bx-ax); 
y= 0:dy:(by-ay);
x = repmat(x,1,1,steptime+3);
y = repmat(y,1,1,steptime+3);     
%% Bottom BC

        for t = 1:steptime+1
            for j = 2:1:stepy+3
               u(1,j,t) = j.*(j-ax).^2; 
            end
        end
        

%% Upper BC
        
        for t = 1:steptime+1
            for j = 2:stepy+3
               u(stepy+3,j,t) = (i-ax)^2 * cos((pi*(j-ax)/(bx-ax)));
            end
        end


%% Right BC
        gabx = bx*(bx-ax)^2;
        fabx = (bx-ax)^2 * cos((pi*(bx-ax)/(bx-ax)));

        for t = 1:steptime+1
            for i = 2:stepx+3
               u(i,stepx+3,t) = gabx + (i-ay)/(by-ay) * (fabx - gabx);
            end
        end
        
if exist( 'checkpoint.mat','file' ) % If a checkpoint file exists, load it
    fprintf('Checkpoint file found - Loading\n');
    load('checkpoint.mat')
end        
%% Explicit Calc

        for t = 1:1:steptime+3
            for i = 2:1:stepy+2
              for j = 2:1:stepx+2
                           
                    u(i,j,t+1) = u(i,j,t) + (dt/dx^2)*(u(i,j+1,t) - 2*u(i,j,t) + u(i,j-1,t)) + (dt/dy^2)*(u(i+1,j,t) - 2*u(i,j,t) + u(i-1,j,t));  
              end
                 
        for r = 2:stepx+1             % Neumann BC
            u(stepy+2,r,t) = u(stepy+1,r,t);
        end  
            end
%Used to plot the solution
%             h = surf(x(:,:,1),y(:,:,1),u(:,:,t),'EdgeColor','none');
%             colorbar
%             shading faceted
%             title({['time = ',num2str(t*dt)]})
%             xlabel('x')
%             ylabel('y')
%             zlabel('Solution ')
%             drawnow; 
%             refreshdata(h)
    if mod(t,t/4)==0
        %save checkpoint   
        fprintf('Saving checkpoint\n');
        save('checkpoint.mat');
    end
        end

        

