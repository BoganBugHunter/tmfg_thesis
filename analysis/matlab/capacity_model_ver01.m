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











%1D area versus position

NozzleLength = 1;
AreaAtThroat = 0.1;
AreaAtEnds = 0.2;

positions = 0:0.001:NozzleLength;

Areas = (4/NozzleLength^2)*(AreaAtEnds-AreaAtThroat).*positions.^2 - (4/NozzleLength)*(AreaAtEnds-AreaAtThroat).*positions + AreaAtEnds;
    


%critical pressure ratio
PressureRatioCritical = ((g+1)/2)^(g/(1-g));

%choked capacity
CapacityChoked = S*AreaAtThroat*sqrt(   ((g+1)/2)^(2/(1-g))  -  ((g+1)/2)^((1+g)/(1-g))   );

%pressure ratios
PressureRatiosArray = PressureRatioCritical:0.0001:1;




%capacity
Capacities = S*AreaAtThroat.*sqrt(PressureRatiosArray.^(2/g).*(1-PressureRatiosArray.^((g-1)/g)));

%mass flow rate
MassFlowRates = (p0/sqrt(T0))*Capacities;









%%%%The 1D CFD

%throat pressure ratio
r= 1.01*PressureRatioCritical;

%single capacity for testing the solver
G = S*AreaAtThroat*sqrt(r^(2/g)*(1-r^((g-1)/g)));

%single mass flow rate for testing the solver
m = (p0/sqrt(T0))*G;

velocity_solutions = zeros(1,length(positions));
pressure_solutions = zeros(1,length(positions));
temperature_solutions = zeros(1,length(positions));
mach_number_solutions = zeros(1,length(positions));



for i = 1:length(positions)
    
    x = positions(i);
    A = (4/NozzleLength^2)*(AreaAtEnds-AreaAtThroat)*x^2 - (4/NozzleLength)*(AreaAtEnds-AreaAtThroat)*x + AreaAtEnds;
    
    %solution coefficients
    alpha = p0^(1-g)*((T0*m*R)/A)^g;
    beta = p0;
    gamma = 1/(R*T0);
    delta = 0.5*(g-1)/g;
    iota = p0/r0;
    
    possible_velocities = 0:0.05:500;
    
    %the sign of the initial difference between the two sides of the
    %equation - left minus right
    sign_init = sign( (alpha/beta)*0^(-g) - ( gamma*iota - (gamma*delta)*0.^2 )^(g/(g-1) ) );
    
    solution_is_found = false;
    
    for j = 1:length(possible_velocities)
        
        v = possible_velocities(j);
       
        %the sign of the current difference between the two sides of the
        %equation - left minus right
        sign_current = sign( (alpha/beta)*v^(-g) - ( gamma*iota - (gamma*delta)*v.^2 )^(g/(g-1) ) );
        
        if sign_current ~= sign_init && solution_is_found == false
            velocity_solutions(i) = v;
            pressure_solutions(i) = p0^(1-g)*((T0*m*R)/(v*A))^g;
            temperature_solutions(i) = v*pressure_solutions(i)*A/(m*R);
            mach_number_solutions(i) = v/sqrt(g*R*temperature_solutions(i));
            solution_is_found = true;
        end
        
    end
    
    
    
end



%equation of the curve of the nozzle in 1D
%Areas = (4/NozzleLength^2)*(AreaAtEnds-AreaAtThroat).*xs.^2 - (4/NozzleLength)*(AreaAtEnds-AreaAtThroat).*xs + AreaAtEnds;


%single form for testing the solver
x = 0.6;
A = (4/NozzleLength^2)*(AreaAtEnds-AreaAtThroat)*x^2 - (4/NozzleLength)*(AreaAtEnds-AreaAtThroat)*x + AreaAtEnds;

%example throat pressure ratio
r = 0.8;

%single capacity for testing the solver
G = S*AreaAtThroat*sqrt(r^(2/g)*(1-r^((g-1)/g)));

%single mass flow rate for testing the solver
m = (p0/sqrt(T0))*G;



%solution coefficients
alpha = p0^(1-g)*((T0*m*R)/A)^g;
beta = p0;
gamma = 1/(R*T0);
delta = 0.5*(g-1)/g;
iota = p0/r0;

%possible values for velocity
vs = 100:0.01:200;

%the two sides of the equation
ys_left = (alpha/beta).*vs.^(-g);
ys_right = ( gamma*iota - (gamma*delta).*vs.^2 ).^(g/(g-1));





%0D capacity plot
figure(1)
hold on
xline(PressureRatioCritical,'--r');
yline(CapacityChoked,'--b');
plot(PressureRatiosArray,Capacities,'g','LineWidth',2);


%1D area vesus position
figure(2)
hold on
plot(positions,Areas,'g','LineWidth',2);

%solution test
figure(3)
hold on
plot(vs,ys_left,'b','LineWidth',2);
plot(vs,ys_right,'r','LineWidth',2);

%plot of position along nozzle vs velocity
figure(4)
hold on
plot(positions,velocity_solutions,'b','LineWidth',2);

%plot of position along nozzle vs velocity
figure(5)
hold on
plot(positions,mach_number_solutions,'k','LineWidth',2);






