close all 
clear variables


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







A_throat = 1; %for now



pressures = 0.1*p0:1:0.9*p0;


areas = zeros(size(pressures));
areas_check = zeros(size(pressures));
Ms = zeros(size(pressures));
numerators = zeros(size(pressures));
denominators = zeros(size(pressures));

for i = 1:size(pressures,2)
    
    p = pressures(i);
    
    constant = sqrt( ((g+1)/2)^((g+1)/(1-g)) * (g-1)/2 );

    numerator = (p0/p)^((g+1)/g);
    numerators(i) = numerator;
    
    denominator = (p0/p)^((g-1)/g) - 1;
    denominators(i) = denominator;

    A = A_throat*constant*sqrt(numerator/denominator);
    areas(i) = A;
    
    
    
    %check all this using the basic components
    
    M = sqrt( (2/(g-1)) * ( (p/p0)^((1-g)/g) -1 ) );
    Ms(i) = M;
    
    Acheck = A_throat * (1/M) * ( (g+1)/2 )^( (g+1)/(2-2*g) ) * ( 1 + ( (g-1)/2 )*M^2 )^( (g+1)/(2*g-2) );
    areas_check(i) = Acheck;
    
    

end



figure(1)
hold on
xline(PressureRatioCritical)
plot(pressures/p0,areas);

figure(2)
hold on
xline(PressureRatioCritical)
plot(pressures/p0,areas_check);

figure(3)
hold on
xline(PressureRatioCritical)
plot(pressures/p0,Ms);









