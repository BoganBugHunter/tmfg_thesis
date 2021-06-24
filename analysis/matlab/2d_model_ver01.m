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

%initial grid sizes
elementSizeAlong = 0.001;
elementSizeAcross = 0.0001;

%number of points in each direction
elementNumberAlong = initialStreamlineLength/elementSizeAlong;
elementNumberAcross = initialPassageWidth/elementSizeAcross;

%bottom left corner
point0 = [0,-initialPassageWidth/2];

%generate the mesh - a column is a streamline
Ax = zeros(elementNumberAlong,elementNumberAcross);
Ay = zeros(elementNumberAlong,elementNumberAcross);
for i = 1:elementNumberAlong
    for j = 1:elementNumberAcross
        Ax(i,j) = point0(1)+(i-1)*elementSizeAlong;
        Ay(i,j) = point0(2)+(j-1)*elementSizeAcross;
    end
    
end

%transform the mesh
Bx = zeros(elementNumberAlong,elementNumberAcross);
By = zeros(elementNumberAlong,elementNumberAcross);
for i = 1:elementNumberAlong
    for j = 1:elementNumberAcross
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
Bangles = zeros(elementNumberAlong,elementNumberAcross);
Bwidths = zeros(elementNumberAlong,elementNumberAcross);
Blengths = zeros(elementNumberAlong,elementNumberAcross);
for i = 1:elementNumberAlong-1
    for j = 1:elementNumberAcross-1
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

%get the local radii of curvature
Bturnings = zeros(elementNumberAlong,elementNumberAcross);
Bradii = zeros(elementNumberAlong,elementNumberAcross);
for i = 1:elementNumberAlong-2
    for j = 1:elementNumberAcross-2
        localTurning = abs( Bangles(i+1,j) - Bangles(i,j) );
        localRadius = Blengths(i,j)/localTurning;
        Bturnings(i,j) = localTurning;
        Bradii(i,j) = localRadius;
    end
end

%trim all the matrices to the size of the smallest ones
Bx = Bx(1:i,1:j);
By = By(1:i,1:j);
Bangles = Bangles(1:i,1:j);
Bwidths = Bwidths(1:i,1:j);
Blengths = Blengths(1:i,1:j);
Bturnings = Bturnings(1:i,1:j);
Bradii = Bradii(1:i,1:j);

%find the choked mass flow rate of each streamline
BthroatAreas = zeros(1,elementNumberAcross);
BmassFlowRates = zeros(1,elementNumberAcross);
for j = 1:elementNumberAcross-2
    streamlineAreas = Bwidths(:,j);
    streamlineThroatArea = min(streamlineAreas);
    BthroatAreas(1,j) = streamlineThroatArea;
    BmassFlowRates(1,j) = (p0/sqrt(T0))*S*streamlineThroatArea*sqrt( ((g+1)/2)^(2/(1-g)) - ((g+1)/2)^((1+g)/(1-g)) );
end
totalThroatArea = sum(BthroatAreas);
totalMassFlowRate = sum(BmassFlowRates);

%calculate velocities at each node
for i = 1:elementNumberAlong-2
    for j = 1:elementNumberAcross-2
        
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



