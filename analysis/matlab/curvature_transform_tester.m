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
rearLoadingFactor = 0.25;
nozzleTurningAngle = 70/360*2*pi;
nozzleRadiusOfCurvature = NozzleLength/nozzleTurningAngle;

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
%first in terms of area changes
for i = 1:numberOfNodesAlong
    for j = 1:numberOfNodesAcross
        area_function = (4/initialStreamlineLength^2)*(AreaAtEnds-AreaAtThroat)...
            *Ax(i,j)^2 - (4/initialStreamlineLength)*(AreaAtEnds-AreaAtThroat)*Ax(i,j) + AreaAtEnds;
        %the transformation matrix for each coordinate
        M = [ 1-Ay(i,j) Ax(i,j) ; Ay(i,j) 0.5*area_function/(initialPassageWidth/2)-Ax(i,j) ];
        
        %MthroatLocation = [ (1-rearLoadingFactor) rearLoadingFactor*initialStreamlineLength/Ay(i,j) ; 0 1 ];
        MthroatLocation = [ 1 + 4*rearLoadingFactor/initialStreamlineLength - Ax(i,j)*4*rearLoadingFactor/initialStreamlineLength^2 0 ; 0 1  ];
        
        M2 = [ (nozzleRadiusOfCurvature/Ax(i,j))*sin(Ax(i,j)*nozzleTurningAngle/NozzleLength)...
            sin(Ax(i,j)*nozzleTurningAngle/NozzleLength) ;...
            (nozzleRadiusOfCurvature/Ax(i,j))*(cos(Ax(i,j)*nozzleTurningAngle/NozzleLength) -1)...
            cos(Ax(i,j)*nozzleTurningAngle/NozzleLength) ];
        
        localVector = [Ax(i,j),Ay(i,j)]';
        outputVector = M2*MthroatLocation*M*localVector;
        Bx(i,j) = outputVector(1);
        By(i,j) = outputVector(2);
    end 
end
%then in terms of a turning angle
%for i = 1:numberOfNodesAlong
%    for j = 1:numberOfNodesAcross
%        %the transformation matrix for each coordinate
%        M2 = [ (nozzleRadiusOfCurvature/Bx(i,j))*sin(Bx(i,j)*nozzleTurningAngle/NozzleLength)...
%            sin(Bx(i,j)*nozzleTurningAngle/NozzleLength) ;...
%            (nozzleRadiusOfCurvature/Bx(i,j))*(cos(Bx(i,j)*nozzleTurningAngle/NozzleLength) -1)...
%            cos(Bx(i,j)*nozzleTurningAngle/NozzleLength) ];
%        localVector = [Bx(i,j),By(i,j)]';
%        outputVector = M2*localVector;
%        Bx(i,j) = outputVector(1);
%        By(i,j) = outputVector(2);
%    end
%end



%plot of the initial square mesh
figure(1)
hold on
plot(Ax,Ay);

%plot of the transformed mesh
figure(2)
hold on
axis equal
plot(Bx,By);