%ver01 goes as far as computing the pressure gradient field as a scalar
%field. It needs to be a vector field, so it can be integrated to give a
%scalar field of pressure corrections due to turning, which can then be
%superposed on the existing scalar pressure field.

%ver02 gets the pressure gradient field as a vector field and integrates it
%to get a scalar pressure field. This requires judicious selection of which
%streamline is to be the baseline (no pressure correction along that
%streamline).

%ver03 includes a second transformation matrix to apply turning to the
%flow, along with an improved method for calculating streamline widths.

%ver04 incorporates another transformation matrix to move the throat back
%and forth.

%ver05 uses more realistic boundary conditions to try to make the Mach 1
%line distortion more clear.



close all;
clear variables;



%fudge factor to exaggerate the M1 line distortion
fudgeFactor = 1; %1 is nominal

%physical constants
R = 287.058; %[Joules/(kg.K)]
cp = 1121; %[Joules/(kg.K)]
cv = 834; %[Joules/(kg.K)]
g = cp/cv; 

%atmospheric conditions
p0 = 4.3256e6; %[Pascals] inlet total pressure from old thesis
T0 = 300; %[Kelvin]

%emergent intitial conditions
r0 = p0/(T0*R); %density

%scale constant
S = sqrt(2*g/(R*(g-1)));

%critical pressure ratio
PressureRatioCritical = ((g+1)/2)^(g/(1-g));



%nozzle shape
NozzleLength = 0.1;
AreaAtThroat = 0.01;
AreaAtEnds = 0.04;
rearLoadingFactor = -0.016;
%nozzleTurningAngle = 76.891/360*2*pi;
nozzleTurningAngle = 120/360*2*pi;
%nozzleTurningAngle = pi;
%nozzleTurningAngle = 0.00001;
nozzleRadiusOfCurvature = NozzleLength/nozzleTurningAngle;

%inital dimensions
initialStreamlineLength = NozzleLength;
initialPassageWidth = AreaAtEnds;

%number of elements in each direction
elementsAlong = 100;
elementsAcross = 100;
middleStreamline = round(elementsAcross/2) + 1;

%number of nodes in each direction
countAlong = elementsAlong+1;
countAcross = elementsAcross+1;

%grid sizes
gridSizeAlong = initialStreamlineLength/elementsAlong;
gridSizeAcross = initialPassageWidth/elementsAcross;



%generate the mesh - a column is a streamline
spawnPoint = [0,-initialPassageWidth/2];
meshX = zeros(countAlong,countAcross);
meshY = zeros(countAlong,countAcross);
for j = 1:countAcross
    for i = 1:countAlong
        meshX(i,j) = spawnPoint(1)+(i-1)*gridSizeAlong;
        meshY(i,j) = spawnPoint(2)+(j-1)*gridSizeAcross;
    end
end


%transform the mesh
x = zeros(size(meshX));
y = zeros(size(meshY));
for i = 1:countAlong
    for j = 1:countAcross
        %the parabolic function describing the area variation along the nozzle
        area_function = (4/initialStreamlineLength^2)*(AreaAtEnds-AreaAtThroat)...
            *meshX(i,j)^2 - (4/initialStreamlineLength)*(AreaAtEnds-AreaAtThroat)*meshX(i,j) + AreaAtEnds;
        
        
        %the transformation matrices
        
        Marea = [ 1-meshY(i,j) meshX(i,j) ; meshY(i,j) 0.5*area_function/(initialPassageWidth/2)-meshX(i,j) ];
        
        MthroatLocation = [ 1 + 4*rearLoadingFactor/initialStreamlineLength... 
            - meshX(i,j)*4*rearLoadingFactor/initialStreamlineLength^2 0 ; 0 1  ];
        
        Mturn = [ (nozzleRadiusOfCurvature/meshX(i,j))*sin(meshX(i,j)*nozzleTurningAngle/NozzleLength)...
            sin(meshX(i,j)*nozzleTurningAngle/NozzleLength) ;...
            (nozzleRadiusOfCurvature/meshX(i,j))*(cos(meshX(i,j)*nozzleTurningAngle/NozzleLength) -1)...
            cos(meshX(i,j)*nozzleTurningAngle/NozzleLength) ];
        %apply the transformations
        localVector = [meshX(i,j),meshY(i,j)]';
        outputVector = Mturn*MthroatLocation*Marea*localVector;
        x(i,j) = outputVector(1);
        y(i,j) = outputVector(2);
    end 
end
%start of streamline cases
for j = 1:countAcross
    x(1,j) = meshX(1,j);
    y(1,j) = meshY(1,j);
end


%find the line of minimum area
possibleDistances = zeros(size(x,1),size(x,1));
for i = 1:countAlong
    point_on_bottom_x = x(i,1);
    point_on_bottom_y = y(i,1);
    for j = 1:countAlong
        point_on_top_x = x(j,size(x,1));
        point_on_top_y = y(j,size(y,1));
        possibleDistances(i,j) = sqrt( (point_on_bottom_x-point_on_top_x)^2 +...
            (point_on_bottom_y-point_on_top_y)^2 );
    end
end
[min1,index1] = min(possibleDistances);
[min2,index2] = min(min1);
minimumDistance = possibleDistances(index1(index2),index2);
minimumAreaLineX = [ x(index1(index2),1) ; x(index1(index2),size(x,1)) ];
minimumAreaLineY = [ y(index1(index2),1) ; y(index1(index2),size(y,1))  ];



%get the x and y components of the flow
widths = zeros(size(x));
xComponents = zeros(size(x));
yComponents = zeros(size(x));
for j = 1:countAcross-1
    for i = 1:countAlong-1
        streamwiseVector = [ x(i+1,j)-x(i,j) , y(i+1,j)-y(i,j) ]';
        spanwiseVector = [ x(i,j+1)-x(i,j) , y(i,j+1)-y(i,j) ]';
        streamlineWidth = sqrt( norm(spanwiseVector)^2 - (dot(streamwiseVector,spanwiseVector))^2 );
        %save them in matrices
        widths(i,j) = streamlineWidth;
        xComponents(i,j) = streamwiseVector(1);
        yComponents(i,j) = streamwiseVector(2);
    end
end
%side wall of passage cases
for i = 1:countAlong
    widths(i,countAcross) = widths(i,countAcross-1);
    xComponents(i,countAcross) = xComponents(i,countAcross-1);
    yComponents(i,countAcross) = yComponents(i,countAcross-1);
end
%end of streamline cases
for j = 1:countAcross
    widths(countAlong,j) = widths(countAlong-1,j);
    xComponents(countAlong,j) = xComponents(countAlong-1,j);
    yComponents(countAlong,j) = yComponents(countAlong-1,j);
end


%get the local radii of curvature as a vector field
radiiXComponents = zeros(countAlong,countAcross);
radiiYComponents = zeros(countAlong,countAcross);
for j = 1:countAcross
    for i = 1:countAlong-2
        radiiXComponents(i,j) = ( xComponents(i,j)*xComponents(i+1,j) + yComponents(i,j)*yComponents(i+1,j) )...
            /( xComponents(i+1,j) - xComponents(i,j)/yComponents(i,j)*yComponents(i+1,j) );
        radiiYComponents(i,j) = -radiiXComponents(i,j)*xComponents(i,j)/yComponents(i,j);
        if isnan(radiiXComponents(i,j))
            radiiXComponents(i,j) = inf;
        end
        if isnan(radiiYComponents(i,j))
            radiiYComponents(i,j) = inf;
        end
    end
end
%end of streamline cases
for j = 1:countAcross
    radiiXComponents(countAlong,j) = radiiXComponents(countAlong-2,j);
    radiiXComponents(countAlong-1,j) = radiiXComponents(countAlong-2,j);
    radiiYComponents(countAlong,j) = radiiYComponents(countAlong-2,j);
    radiiYComponents(countAlong-1,j) = radiiYComponents(countAlong-2,j);
end


%find the choked mass flow rate of each streamline
throatAreas = zeros(1,countAcross);
massFlowRates = zeros(1,countAcross);
for j = 1:countAcross
    streamlineAreas = widths(:,j);
    throatAreas(1,j) = min(streamlineAreas);
    massFlowRates(1,j) = (p0/sqrt(T0))*S*throatAreas(1,j)*sqrt( ((g+1)/2)^(2/(1-g)) - ((g+1)/2)^((1+g)/(1-g)) );
end
totalThroatArea = sum(throatAreas);
totalMassFlowRate = sum(massFlowRates);


%calculate the speeds at each node using the solver function
speeds = zeros(size(x,1),size(x,2));
pressures = zeros(size(x,1),size(x,2));
temperatures = zeros(size(x,1),size(x,2));
densities = zeros(size(x,1),size(x,2));
machNumbers = zeros(size(x,1),size(x,2));
speedsOfSound = zeros(size(x,1),size(x,2));
for j = 1:countAcross
    [speeds(:,j),pressures(:,j),temperatures(:,j),densities(:,j),machNumbers(:,j),speedsOfSound(:,j)]...
        = streamlineSolverVer02(widths(:,j),massFlowRates(1,j));
end


%compute the magnitude of the local pressure gradient due to turning at each node
pressureGradientsX = zeros(size(x,1),size(x,2));
pressureGradientsY = zeros(size(x,1),size(x,2));
for j = 1:countAcross
    for i = 1:countAlong
        %get the scalar values for velocity and radius of curvature
        velocityMagnitude = speeds(i,j);
        radiusMagnitude = sqrt( radiiXComponents(i,j)^2 + radiiYComponents(i,j)^2 );
        %pressure gradient due to turning is always in the direction of the radius vector
        pressureGradientX = densities(i,j)*((velocityMagnitude/radiusMagnitude)^2)*radiiXComponents(i,j);
        pressureGradientY = densities(i,j)*((velocityMagnitude/radiusMagnitude)^2)*radiiYComponents(i,j);
        if radiusMagnitude == inf
            pressureGradientX = 0;
            pressureGradientY = 0;
        end
        pressureGradientsX(i,j) = pressureGradientX;
        pressureGradientsY(i,j) = pressureGradientY;
    end
end


%integrate the pressure gradient vector field along parallel vertical lines, to give
%the pressure correction scalar field
integrands = zeros(size(x,1),size(x,2));
pressureCorrections = zeros(size(x,1),size(x,2));
for i = 1:countAlong-1
    pressureCorrection = 0;
    pressureCorrections(i,1) = 0; 
    for j = 2:countAcross-1
        verticalNodeSeparation = yComponents(i,j+1) - yComponents(i,j);
        integrand = pressureGradientsY(i,j)*verticalNodeSeparation;
        integrands(i,j) = integrand;
        pressureCorrection = pressureCorrection + integrand;
        pressureCorrections(i,j) = fudgeFactor*pressureCorrection;
    end
end
%end of streamline cases
for j = 1:countAcross
    integrands(countAlong,j) = integrands(countAlong-1,j);
    pressureCorrections(countAlong,j) = pressureCorrections(countAlong-1,j);
end
%side wall cases
for i = 1:countAlong
    integrands(i,countAcross-1) = integrands(i,countAcross-2);
    integrands(i,countAcross) = integrands(i,countAcross-1);
    pressureCorrections(i,countAcross) = pressureCorrections(i,countAcross-1);
end
%normalise against the DC offset from a suitable streamline (ie get the constant
%of integration
pressureCorrectionsNormalised = zeros(size(x,1),size(x,2));
for j = 1:countAcross
    for i = 1:countAlong
        %test = pressureCorrections(i,middleStreamline); %should it always be the middle streamline?
        test = pressureCorrections(i,size(x,1));
        %test = pressureCorrections(i,1);
        pressureCorrectionsNormalised(i,j) = pressureCorrections(i,j) - test;
    end
end
%apply the pressure corrections, remembering what points they are
%referenced to
pressuresRevised = zeros(size(x,1),size(x,2));
for j = 1:countAcross
    for i = 1:countAlong
        pressuresRevised(i,j) = pressures(i,j) + pressureCorrectionsNormalised(i,j);
    end
end


%turn the new pressures into new Mach numbers
machNumbers2 = zeros(size(x,1),size(x,2));
for j = 1:countAcross
    for i = 1:countAlong
        machNumbers2(i,j) = sqrt( (2/(g-1)) * ( (pressuresRevised(i,j)/p0)^((1-g)/g) - 1 ) );
    end
end

%turn the new Mach numbers into new areas
areasRevised = zeros(size(x,1),size(x,2));
for j = 1:countAcross
   %get this streamline's current throat area
   streamlineAreas = widths(:,j);
   throatArea = min(streamlineAreas);
   %calculate the area ratio of every point in the field and multiply by
   %the throat area
   for i = 1:countAlong
       areasRevised(i,j) = (throatArea/machNumbers2(i,j)) * ( (g+1)/2 /...
           ( 1 + (g-1)*machNumbers2(i,j)^2/2 ) )^((g+1)/2-2*g);
   end
end

%use the new areas to re-run the solver function
for j = 1:countAcross
    streamlineAreas = areasRevised(:,j);
    throatAreas(1,j) = min(streamlineAreas);
    massFlowRates(1,j) = (p0/sqrt(T0))*S*throatAreas(1,j)*sqrt( ((g+1)/2)^(2/(1-g)) - ((g+1)/2)^((1+g)/(1-g)) );
end
totalThroatArea2 = sum(throatAreas);
totalMassFlowRate2 = sum(massFlowRates);

for j = 1:countAcross
    [speeds(:,j),pressures(:,j),temperatures(:,j),densities(:,j),machNumbers(:,j),speedsOfSound(:,j)]...
        = streamlineSolverVer02(widths(:,j),massFlowRates(1,j));
end









%turn the new Mach numbers into new areas - measure the new capacity by measuring
%the areas where M = 1



%turn the pressure corrections into area corrections

%apply the area corrections

%if I transform all the areas in the flow field, it is likely that the
%smallest areas will remain the smallest areas. But there may be some
%places where the location of choking is different because the pressure
%correction is sufficient to make some other area transform into the new
%smallest area

%run the new areas and mass flow rates through the streamline solver to get
%the new quantities along each streamline

%use this to get a capacity solution

%investigate the worth of looping round this approach to convergence?







%plot of the initial square mesh
figure(1)
hold on
axis equal
plot(meshX,meshY);

%plot of the transformed mesh
figure(2)
hold on
axis equal
plot(x,y);

%contours of local area
figure(3)
hold on
axis equal
contourf(x,y,widths,20, 'LineColor', 'none');

%contours of local angle
%figure(4)
%hold on
%contourf(Bx,By,Bangles,20, 'LineColor', 'none');

%contours of local turning
%figure(5)
%hold on
%contourf(Bx,By,Bturnings,20, 'LineColor', 'k');

%contours of local radius of curvature
%figure(6)
%hold on
%contourf(Bx,By,Bradii,20, 'LineColor', 'none');

%plot of position along nozzle vs velocity
figure(7)
hold on
title('Flow velocity vs position')
ylabel('Velocity, m/s')
xlabel('Position')
set(gca,'FontSize',14)
plot(speeds(:,middleStreamline),'k','LineWidth',2);

%contours of velocity
figure(8)
hold on
axis equal
contourf(x,y,speeds,40, 'LineColor', 'k', 'ShowText','on');

%contours of pressure gradient
figure(9)
hold on
axis equal
contourf(x,y,sqrt(pressureGradientsX^2+pressureGradientsY^2),100, 'LineColor', 'none');

%contours of pressure correction
figure(12)
hold on
axis equal
contourf(x,y,pressureCorrectionsNormalised,100, 'LineColor', 'none');

%contours of revised pressure
figure(13)
hold on
axis equal
contourf(x,y,pressuresRevised,100, 'LineColor', 'none');

%contours of Mach number
figure(14)
hold on
axis equal
contourf(x,y,machNumbers,200, 'LineColor', 'none', 'ShowText','off');
contour(x,y,machNumbers,[0,0.5,1.0,1.5], 'LineColor', 'k', 'ShowText','on');
plot(minimumAreaLineX,minimumAreaLineY,'m')





