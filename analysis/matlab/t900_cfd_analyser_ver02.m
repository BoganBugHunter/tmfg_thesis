close all;
clear variables;

load workspaces/T900_2021.mat;

% the spacing in metres between each NGV at mid-span radius
periodic_spacing = 0.05881; 

%load data from the imported workspacew
x = KSZ03_xvel_yvel_mach(:,2);
y = KSZ03_xvel_yvel_mach(:,3);
xvel = KSZ03_xvel_yvel_mach(:,4);
yvel = KSZ03_xvel_yvel_mach(:,5);
M = KSZ03_xvel_yvel_mach(:,6);

%make it 2 vanes instead of 1
x_repeated = x;
y_repeated = y;
for n = 1:length(y)
   y_repeated(n) = y(n) + periodic_spacing;
end
x = cat(1,x,x_repeated);
y = cat(1,y,y_repeated);
M = cat(1,M,M);
xvel = cat(1,xvel,xvel);
yvel = cat(1,yvel,yvel);

%create a uniform quadrilateral mesh onto which the existing data will be
%interpolated
[x_grid,y_grid] = meshgrid(-0.03:0.0001:-0.009,0.050:0.0001:0.070);
mach_grid = griddata(x, y, M, x_grid, y_grid);
xvel_grid = griddata(x, y, xvel, x_grid, y_grid);
yvel_grid = griddata(x, y, yvel, x_grid, y_grid);
mach_grid_size = size(mach_grid);
height = mach_grid_size(1);
width = mach_grid_size(2);

%replace NaN values with 0
for i = 1:height
    for j = 1:width
        if isnan(mach_grid(i,j))
           mach_grid(i,j) = 0; 
        end
    end
end

%find the M=1 line
M1_positions_and_velocities = zeros(1,4);
for j = 1:width
    %go along a vertical strip until M > 1
    for i = 2:height
        if (mach_grid(i-1,j) > 1 && mach_grid(i,j) < 1) ||...
                (mach_grid(i-1,j) < 1 && mach_grid(i,j) > 1)
            M1_position_and_velocity = [x_grid(i,j),y_grid(i,j),...
                xvel_grid(i,j),yvel_grid(i,j)];
            M1_positions_and_velocities =...
                cat(1,M1_positions_and_velocities,...
                M1_position_and_velocity);
        end
    end
end
M1_positions_and_velocities = M1_positions_and_velocities(2:end,:);

%add two columns to M1_positions which is are x and y coordinates tilted at
%an angle approximately orthogonal to the sonic line so they can be sorted 
%along this direction
theta = 2*pi*-30/360;
R = [cos(theta) -sin(theta) ; sin(theta) cos(theta)];
M1_positions_and_velocities = cat(2,M1_positions_and_velocities,...
    zeros(size(M1_positions_and_velocities,1),2));
for i = 1:size(M1_positions_and_velocities,1)
    outputVector = R*[M1_positions_and_velocities(i,1),...
       M1_positions_and_velocities(i,2)]';
    M1_positions_and_velocities(i,5) = outputVector(1);
    M1_positions_and_velocities(i,6) = outputVector(2);
end
M1_positions_and_velocities = sortrows(M1_positions_and_velocities,5);

%trim the ends so spurious bits aren't included
M1_positions_and_velocities = M1_positions_and_velocities(1:end,:);

%integrate along the M1 line
theta = 2*pi*90/360;
R = [cos(theta) -sin(theta) ; sin(theta) cos(theta)];
integral = 0;
for i = 2:size(M1_positions_and_velocities,1)
    localVector = [M1_positions_and_velocities(i,1),...
        M1_positions_and_velocities(i,2)]-...
        [M1_positions_and_velocities(i-1,1),...
        M1_positions_and_velocities(i-1,2)];
    localUnitVector = localVector/norm(localVector);
    localOrthogonalUnitVector = R*localUnitVector';
    %find the velocity vector local to this point on the line
    localVelocity = [M1_positions_and_velocities(i,3),...
        M1_positions_and_velocities(i,4)];
    localVelocityUnitVector = localVelocity/norm(localVelocity);
    integrand = dot(localOrthogonalUnitVector,...
        localVelocityUnitVector)*norm(localVector);
    if isnan(integrand)
    else
    integral = integral + integrand;
    end
end

%plot the sonic line shape over the mesh nodes
figure(1)
plot(x,y,'k.')
hold on
plot(M1_positions_and_velocities(:,1),...
    M1_positions_and_velocities(:,2),'r-')
axis equal
xlim([-0.03 -0.009])
ylim([0.050 0.070])