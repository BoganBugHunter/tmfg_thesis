%ver01 goes as far as computing the pressure gradient field as a scalar
%field. It needs to be a vector field, so it can be integrated to give a
%scalar field of pressure corrections due to turning, which can then be
%superposed on the existing scalar pressure field.



close all;
clear variables;



%physical constants
R = 287.058; %[Joules/(kg.K)]
cp = 1121; %[Joules/(kg.K)]
cv = 834; %[Joules/(kg.K)]
g = cp/cv; 

%atmospheric conditions
r0 = 1.225; %[kg/m3]
p0 = 101325; %[Pascals]

%emergent intitial conditions
T0 = p0/(r0*R);

%scale constant
S = sqrt(2*g/(R*(g-1)));

%critical pressure ratio
PressureRatioCritical = ((g+1)/2)^(g/(1-g));



%nozzle shape
NozzleLength = 1;
AreaAtThroat = 0.01;
AreaAtEnds = 0.1;

%inital dimensions
initialStreamlineLength = NozzleLength;
initialPassageWidth = AreaAtEnds;

%number of elements in each direction
numberOfElementsAlong = 200;
numberOfElementsAcross = 100;

%number of nodes in each direction
numberOfNodesAlong = numberOfElementsAlong+1;
numberOfNodesAcross = numberOfElementsAcross+1;

%grid sizes
elementSizeAlong = initialStreamlineLength/numberOfElementsAlong;
elementSizeAcross = initialPassageWidth/numberOfElementsAcross;



%generate the mesh - a column is a streamline
point0 = [0,-initialPassageWidth/2];
Ax = zeros(numberOfNodesAlong,numberOfNodesAcross);
Ay = zeros(numberOfNodesAlong,numberOfNodesAcross);
for j = 1:numberOfNodesAcross
    for i = 1:numberOfNodesAlong
        Ax(i,j) = point0(1)+(i-1)*elementSizeAlong;
        Ay(i,j) = point0(2)+(j-1)*elementSizeAcross;
    end
end



%transform the mesh
Bx = zeros(size(Ax));
By = zeros(size(Ay));
for i = 1:numberOfNodesAlong
    for j = 1:numberOfNodesAcross
        area_function = (4/initialStreamlineLength^2)*(AreaAtEnds-AreaAtThroat)*Ax(i,j)^2 - (4/initialStreamlineLength)*(AreaAtEnds-AreaAtThroat)*Ax(i,j) + AreaAtEnds;
        %the transformation matrix for each coordinate
        M = [ 1-Ay(i,j) Ax(i,j) ; Ay(i,j) 0.5*area_function/(initialPassageWidth/2)-Ax(i,j) ];
        localVector = [Ax(i,j),Ay(i,j)]';
        outputVector = M*localVector;
        Bx(i,j) = outputVector(1);
        By(i,j) = outputVector(2);
    end
    
end


%get the matrices of turning angles and streamline widths
Bangles = zeros(size(Bx));
Bwidths = zeros(size(Bx));
Blengths = zeros(size(Bx));
for j = 1:numberOfNodesAcross-1
    for i = 1:numberOfNodesAlong-1
        %calculate each length and angle
        horizontalTravel = Bx(i+1,j) - Bx(i,j);
        verticalTravel = By(i+1,j) - By(i,j);
        verticalSeparation = By(i,j+1) - By(i,j);
        segmentLength = sqrt(horizontalTravel^2 + verticalTravel^2);
        angleFromVertical = acos(verticalTravel/segmentLength);
        streamlineWidth = verticalSeparation*sin(angleFromVertical);
        %save them in matrices
        Bangles(i,j) = angleFromVertical;
        Bwidths(i,j) = streamlineWidth;
        Blengths(i,j) = segmentLength;
    end
end
%side wall of passage cases
for i = 1:numberOfNodesAlong
    Bangles(i,numberOfNodesAcross) = Bangles(i,numberOfNodesAcross-1);
    Bwidths(i,numberOfNodesAcross) = Bwidths(i,numberOfNodesAcross-1);
    Blengths(i,numberOfNodesAcross) = Blengths(i,numberOfNodesAcross-1);
end
%end of streamline cases
for j = 1:numberOfNodesAcross
    Bangles(numberOfNodesAlong,j) = Bangles(numberOfNodesAlong-1,j);
    Bwidths(numberOfNodesAlong,j) = Bwidths(numberOfNodesAlong-1,j);
    Blengths(numberOfNodesAlong,j) = Blengths(numberOfNodesAlong-1,j);
end


%get the local radii of curvature
Bturnings = zeros(numberOfNodesAlong,numberOfNodesAcross);
Bradii = zeros(numberOfNodesAlong,numberOfNodesAcross);
for j = 1:numberOfNodesAcross
    for i = 1:numberOfNodesAlong-2
        localTurning = abs( Bangles(i+1,j) - Bangles(i,j) );
        localRadius = Blengths(i,j)/localTurning;
        Bturnings(i,j) = localTurning;
        Bradii(i,j) = localRadius;
    end
end
%end of streamline cases
for j = 1:numberOfNodesAcross
    Bturnings(numberOfNodesAlong,j) = Bturnings(numberOfNodesAlong-2,j);
    Bturnings(numberOfNodesAlong-1,j) = Bturnings(numberOfNodesAlong-2,j);
    Bradii(numberOfNodesAlong,j) = Bradii(numberOfNodesAlong-2,j);
    Bradii(numberOfNodesAlong-1,j) = Bradii(numberOfNodesAlong-2,j);
end


%find the choked mass flow rate of each streamline
BthroatAreas = zeros(1,numberOfNodesAcross);
BmassFlowRates = zeros(1,numberOfNodesAcross);
for j = 1:numberOfNodesAcross
    streamlineAreas = Bwidths(:,j);
    streamlineThroatArea = min(streamlineAreas);
    BthroatAreas(1,j) = streamlineThroatArea;
    BmassFlowRates(1,j) = (p0/sqrt(T0))*S*streamlineThroatArea*sqrt( ((g+1)/2)^(2/(1-g)) - ((g+1)/2)^((1+g)/(1-g)) );
end
totalThroatArea = sum(BthroatAreas);
totalMassFlowRate = sum(BmassFlowRates);


%calculate the velocities at each node
Bvelocities = zeros(size(Bx,1),size(Bx,2));
Bpressures = zeros(size(Bx,1),size(Bx,2));
Btemperatures = zeros(size(Bx,1),size(Bx,2));
Bdensities = zeros(size(Bx,1),size(Bx,2));
BmachNumbers = zeros(size(Bx,1),size(Bx,2));
BspeedsOfSound = zeros(size(Bx,1),size(Bx,2));
BareaRatios = zeros(size(Bx,1),size(Bx,2));
for j = 1:numberOfNodesAcross
    streamlineAreas = Bwidths(:,j);
    streamlineMassFlowRate = BmassFlowRates(1,j);
    [streamlineVelocities,streamlinePressures,streamlineTemperatures,streamlineDensities,...
        streamlineMachNumbers,streamlineSpeedsOfSound,streamlineAreaRatios]...
        = streamlineSolver_ver01(streamlineAreas,streamlineMassFlowRate);
    Bvelocities(:,j) = streamlineVelocities;
    Bpressures(:,j) = streamlinePressures;
    Btemperatures(:,j) = streamlineTemperatures;
    Bdensities(:,j) = streamlineDensities;
    BmachNumbers(:,j) = streamlineMachNumbers;
    BspeedsOfSound(:,j) = streamlineSpeedsOfSound;
    BareaRatios(:,j) = streamlineAreaRatios;
end


%calculate the local pressure gradient due to turning at each node
BpressureGradients = zeros(size(Bx,1),size(Bx,2));
for j = 1:numberOfNodesAcross
    for i = 1:numberOfNodesAlong
        BpressureGradients(i,j) = Bdensities(i,j)*Bvelocities(i,j)^2/Bradii(i,j);
    end
end



%plot of the initial square mesh
figure(1)
hold on
plot(Ax,Ay);

%plot of the transformed mesh
figure(2)
hold on
plot(Bx,By);

%contours of local area
figure(3)
hold on
contourf(Bx,By,Bwidths,20, 'LineColor', 'none');

%contours of local angle
figure(4)
hold on
contourf(Bx,By,Bangles,20, 'LineColor', 'none');

%contours of local turning
figure(5)
hold on
contourf(Bx,By,Bturnings,20, 'LineColor', 'k');

%contours of local radius of curvature
figure(6)
hold on
contourf(Bx,By,Bradii,20, 'LineColor', 'none');

%plot of position along nozzle vs velocity
figure(7)
hold on
title('Flow velocity vs position')
ylabel('Velocity, m/s')
xlabel('Position')
set(gca,'FontSize',14)
plot(streamlineVelocities,'k','LineWidth',2);

%contours of velocity
figure(8)
hold on
contourf(Bx,By,Bvelocities,40, 'LineColor', 'k');

%contours of pressure gradient
figure(9)
hold on
contourf(Bx,By,BpressureGradients,100, 'LineColor', 'none');

%contours of density
figure(10)
hold on
contourf(Bx,By,Bdensities,100, 'LineColor', 'none');



