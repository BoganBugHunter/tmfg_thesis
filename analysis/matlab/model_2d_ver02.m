%ver01 goes as far as computing the pressure gradient field as a scalar
%field. It needs to be a vector field, so it can be integrated to give a
%scalar field of pressure corrections due to turning, which can then be
%superposed on the existing scalar pressure field.

%ver02 gets the pressure gradient field as a vector field and integrates it
%to get a scalar pressure field. This requires judicious selection of which
%streamline is to be the baseline (no pressure correction along that
%streamline).



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
        area_function = (4/initialStreamlineLength^2)*(AreaAtEnds-AreaAtThroat)...
            *Ax(i,j)^2 - (4/initialStreamlineLength)*(AreaAtEnds-AreaAtThroat)*Ax(i,j) + AreaAtEnds;
        %the transformation matrix for each coordinate
        M = [ 1-Ay(i,j) Ax(i,j) ; Ay(i,j) 0.5*area_function/(initialPassageWidth/2)-Ax(i,j) ];
        localVector = [Ax(i,j),Ay(i,j)]';
        outputVector = M*localVector;
        Bx(i,j) = outputVector(1);
        By(i,j) = outputVector(2);
    end
    
end


%get the matrices of x and y components of the flow
Bwidths = zeros(size(Bx));
BxComponents = zeros(size(Bx));
ByComponents = zeros(size(Bx));
for j = 1:numberOfNodesAcross-1
    for i = 1:numberOfNodesAlong-1
        %calculate each component, length, and angle
        horizontalTravel = Bx(i+1,j) - Bx(i,j);
        verticalNodeSeparation = By(i+1,j) - By(i,j);
        verticalSeparation = By(i,j+1) - By(i,j);
        segmentLength = sqrt(horizontalTravel^2 + verticalNodeSeparation^2);
        angleFromVertical = acos(verticalNodeSeparation/segmentLength);
        streamlineWidth = verticalSeparation*sin(angleFromVertical);
        %save them in matrices
        Bwidths(i,j) = streamlineWidth;
        BxComponents(i,j) = horizontalTravel;
        ByComponents(i,j) = verticalNodeSeparation;
    end
end
%side wall of passage cases
for i = 1:numberOfNodesAlong
    Bwidths(i,numberOfNodesAcross) = Bwidths(i,numberOfNodesAcross-1);
    BxComponents(i,numberOfNodesAcross) = BxComponents(i,numberOfNodesAcross-1);
    ByComponents(i,numberOfNodesAcross) = ByComponents(i,numberOfNodesAcross-1);
end
%end of streamline cases
for j = 1:numberOfNodesAcross
    Bwidths(numberOfNodesAlong,j) = Bwidths(numberOfNodesAlong-1,j);
    BxComponents(numberOfNodesAlong,j) = BxComponents(numberOfNodesAlong-1,j);
    ByComponents(numberOfNodesAlong,j) = ByComponents(numberOfNodesAlong-1,j);
end


%get the local radii of curvature as a vector field
BradiiXComponents = zeros(numberOfNodesAlong,numberOfNodesAcross);
BradiiYComponents = zeros(numberOfNodesAlong,numberOfNodesAcross);
for j = 1:numberOfNodesAcross
    for i = 1:numberOfNodesAlong-2
        BradiiXComponents(i,j) = ( BxComponents(i,j)*BxComponents(i+1,j) + ByComponents(i,j)*ByComponents(i+1,j) )...
            /( BxComponents(i+1,j) - BxComponents(i,j)/ByComponents(i,j)*ByComponents(i+1,j) );
        BradiiYComponents(i,j) = -BradiiXComponents(i,j)*BxComponents(i,j)/ByComponents(i,j);
        if isnan(BradiiXComponents(i,j))
            BradiiXComponents(i,j) = inf;
        end
        if isnan(BradiiYComponents(i,j))
            BradiiYComponents(i,j) = inf;
        end
    end
end
%end of streamline cases
for j = 1:numberOfNodesAcross
    BradiiXComponents(numberOfNodesAlong,j) = BradiiXComponents(numberOfNodesAlong-2,j);
    BradiiXComponents(numberOfNodesAlong-1,j) = BradiiXComponents(numberOfNodesAlong-2,j);
    BradiiYComponents(numberOfNodesAlong,j) = BradiiYComponents(numberOfNodesAlong-2,j);
    BradiiYComponents(numberOfNodesAlong-1,j) = BradiiYComponents(numberOfNodesAlong-2,j);
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


%calculate the speeds at each node
Bspeeds = zeros(size(Bx,1),size(Bx,2));
Bpressures = zeros(size(Bx,1),size(Bx,2));
Btemperatures = zeros(size(Bx,1),size(Bx,2));
Bdensities = zeros(size(Bx,1),size(Bx,2));
BmachNumbers = zeros(size(Bx,1),size(Bx,2));
BspeedsOfSound = zeros(size(Bx,1),size(Bx,2));
BareaRatios = zeros(size(Bx,1),size(Bx,2));
for j = 1:numberOfNodesAcross
    streamlineAreas = Bwidths(:,j);
    streamlineMassFlowRate = BmassFlowRates(1,j);
    [streamlineSpeeds,streamlinePressures,streamlineTemperatures,streamlineDensities,...
        streamlineMachNumbers,streamlineSpeedsOfSound,streamlineAreaRatios]...
        = streamlineSolver_ver01(streamlineAreas,streamlineMassFlowRate);
    Bspeeds(:,j) = streamlineSpeeds;
    Bpressures(:,j) = streamlinePressures;
    Btemperatures(:,j) = streamlineTemperatures;
    Bdensities(:,j) = streamlineDensities;
    BmachNumbers(:,j) = streamlineMachNumbers;
    BspeedsOfSound(:,j) = streamlineSpeedsOfSound;
    BareaRatios(:,j) = streamlineAreaRatios;
end


%compute the magnitude of the local pressure gradient due to turning at each node
BpressureGradientsScalars = zeros(size(Bx,1),size(Bx,2));
BpressureGradientsX = zeros(size(Bx,1),size(Bx,2));
BpressureGradientsY = zeros(size(Bx,1),size(Bx,2));
for j = 1:numberOfNodesAcross
    for i = 1:numberOfNodesAlong
        %get the scalar values for velocity and radius of curvature
        velocityMagnitude = Bspeeds(i,j);
        radiusMagnitude = sqrt( BradiiXComponents(i,j)^2 + BradiiYComponents(i,j)^2 );
        %pressure gradient is always in the direction of the radius vector
        pressureGradientX = Bdensities(i,j)*((velocityMagnitude/radiusMagnitude)^2)*BradiiXComponents(i,j);
        pressureGradientY = Bdensities(i,j)*((velocityMagnitude/radiusMagnitude)^2)*BradiiYComponents(i,j);
        if radiusMagnitude == inf
            pressureGradientX = 0;
            pressureGradientY = 0;
        end
        BpressureGradientsScalars(i,j) = sqrt( pressureGradientX^2 + pressureGradientY^2 );
        BpressureGradientsX(i,j) = pressureGradientX;
        BpressureGradientsY(i,j) = pressureGradientY;
    end
end


%integrate the pressure gradient vector field along parallel vertical lines, to give
%the pressure correction scalar field
integrands = zeros(size(Bx,1),size(Bx,2));
BpressureCorrections = zeros(size(Bx,1),size(Bx,2));
for i = 1:numberOfNodesAlong-1
    pressureCorrection = 0;
    BpressureCorrections(i,1) = 0; 
    for j = 2:numberOfNodesAcross-1
        verticalNodeSeparation = ByComponents(i,j+1) - ByComponents(i,j);
        integrand = BpressureGradientsY(i,j)*verticalNodeSeparation;
        integrands(i,j) = integrand;
        pressureCorrection = pressureCorrection + integrand;
        BpressureCorrections(i,j) = pressureCorrection;
    end
end
%end of streamline cases
for j = 1:numberOfNodesAcross
    integrands(numberOfNodesAlong,j) = integrands(numberOfNodesAlong-1,j);
    BpressureCorrections(numberOfNodesAlong,j) = BpressureCorrections(numberOfNodesAlong-1,j);
end
%side wall cases
for i = 1:numberOfNodesAlong
    integrands(i,numberOfNodesAcross-1) = integrands(i,numberOfNodesAcross-2);
    integrands(i,numberOfNodesAcross) = integrands(i,numberOfNodesAcross-1);
    BpressureCorrections(i,numberOfNodesAcross) = BpressureCorrections(i,numberOfNodesAcross-1);
end
%normalise against the DC offset from a suitable point (ie get the constant
%of integration
BpressureCorrectionsNormalised = zeros(size(Bx,1),size(Bx,2));
for j = 1:numberOfNodesAcross
    for i = 1:numberOfNodesAlong
        test = BpressureCorrections(i,51);
        BpressureCorrectionsNormalised(i,j) = BpressureCorrections(i,j) - test;
    end
end

%apply the pressure corrections, remembering what points they are
%referenced to
BpressuresRevised = zeros(size(Bx,1),size(Bx,2));
for j = 1:numberOfNodesAcross
    for i = 1:numberOfNodesAlong
        BpressuresRevised(i,j) = Bpressures(i,j) + BpressureCorrections(i,j);
    end
end




%check the radii of curvature are always perpendicular to the streamlines
%the resulting matrix will be full of zeros if they are really
%perpendicular
Btests = zeros(size(Bx,1),size(Bx,2));
for j = 1:numberOfNodesAcross
    for i = 1:numberOfNodesAlong
        Btests(i,j) = BxComponents(i,j)*BpressureGradientsX(i,j) + ByComponents(i,j)*BpressureGradientsY(i,j);
    end
end





%plot of the initial square mesh
%figure(1)
%hold on
%plot(Ax,Ay);

%plot of the transformed mesh
%figure(2)
%hold on
%plot(Bx,By);

%contours of local area
%figure(3)
%hold on
%contourf(Bx,By,Bwidths,20, 'LineColor', 'none');

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
%figure(7)
%hold on
%title('Flow velocity vs position')
%ylabel('Velocity, m/s')
%xlabel('Position')
%set(gca,'FontSize',14)
%plot(streamlineSpeeds,'k','LineWidth',2);

%contours of velocity
%figure(8)
%hold on
%contourf(Bx,By,Bspeeds,40, 'LineColor', 'k');

%contours of pressure gradient
%figure(9)
%hold on
%contourf(Bx,By,BpressureGradientsScalars,100, 'LineColor', 'none');

%contours of the integrand
%figure(10)
%hold on
%contourf(Bx,By,integrands,100, 'LineColor', 'none');

%contours of the integrand
%figure(11)
%hold on
%contourf(Bx,By,BpressureCorrections,100, 'LineColor', 'none');

%contours of pressure correction
figure(12)
hold on
contourf(Bx,By,BpressureCorrectionsNormalised,100, 'LineColor', 'none');

%contours of revised pressure
%figure(13)
%hold on
%contourf(Bx,By,BpressuresRevised,100, 'LineColor', 'none');

%contours of original pressure
%figure(14)
%hold on
%contourf(Bx,By,Bpressures,100, 'LineColor', 'none');





