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

%nozzle shape
NozzleLength = 1;
AreaAtThroat = 0.01;
AreaAtEnds = 0.1;



%%%%Capacity vs pressure ratio for a simple nozzle


%critical pressure ratio
PressureRatioCritical = ((g+1)/2)^(g/(1-g));

%choked capacity
CapacityChoked = S*AreaAtThroat*sqrt(   ((g+1)/2)^(2/(1-g))  -  ((g+1)/2)^((1+g)/(1-g))   );

%pressure ratios
PressureRatiosArray = PressureRatioCritical:0.001:1;

%capacity
Capacities = S*AreaAtThroat.*sqrt(PressureRatiosArray.^(2/g).*(1-PressureRatiosArray.^((g-1)/g)));

%mass flow rate
MassFlowRates = (p0/sqrt(T0))*Capacities;



%%%%The 1D CFD

 
%1D area versus position
positions = 0:0.01:NozzleLength*1;
Areas = (4/NozzleLength^2)*(AreaAtEnds-AreaAtThroat).*positions.^2 - (4/NozzleLength)*(AreaAtEnds-AreaAtThroat).*positions + AreaAtEnds;

%throat pressure ratio
r = PressureRatioCritical;

%mass flow rate
m = (p0/sqrt(T0))*S*AreaAtThroat*sqrt(r^(2/g)*(1-r^((g-1)/g)));

velocity_solutions = zeros(2,length(positions));
pressure_solutions = zeros(1,length(positions));
density_solutions = zeros(1,length(positions));
temperature_solutions = zeros(1,length(positions));
mach_number_solutions = zeros(1,length(positions));
speed_of_sound_solutions = zeros(1,length(positions));
area_ratio_solutions = zeros(1,length(positions));

critical = false;


%maximum speed of the nozzle
starting_guess_high = sqrt(  2*p0/r0*(g/(g-1))  ); 

%minimum speed of the nozzle
starting_guess_low = 0;

%track the total number of iterations
total_iterations = 0;


for i = 1:length(positions)
    
    
    x = positions(i);
    A = (4/NozzleLength^2)*(AreaAtEnds-AreaAtThroat)*x^2 - (4/NozzleLength)*(AreaAtEnds-AreaAtThroat)*x + AreaAtEnds;
    
    
    iterated_solutions = zeros(2,1e6);
    
    
    %start at the starting guesses, but after the first answers, use each
    %velocity to guess the next one downstream
    if i == 1
        iterated_solutions(1,1) = starting_guess_high;
        iterated_solutions(2,1) = starting_guess_low;
    else
        iterated_solutions(1,1) = velocity_solutions(1,i-1);
        iterated_solutions(2,1) = velocity_solutions(2,i-1);
    end
    
    
    residual = starting_guess_high;
    
    n = 2;
    
    while residual > 1e-9

        %this equation is real between m/(r0*A) and infinity
        iterated_solutions(1,n) = sqrt(  2*(p0/r0)*(g/(g-1))*( 1 - (m/(r0*A*iterated_solutions(1,n-1)))^(g-1) )  );
    
        %this equation is real between 0 and the maximum nozzle speed
        iterated_solutions(2,n) = (m/(r0*A))*(  1 - 0.5*(r0/p0)*((g-1)/g)*iterated_solutions(2,n-1)^2  )^(1/(1-g)) ;

        residual = abs( iterated_solutions(2,n) - iterated_solutions(2,n-1) );
    
        solution_high = iterated_solutions(1,n);
        solution_low = iterated_solutions(2,n);
    
        n = n + 1;
        
    end
    
    iterations_to_converge = n-1;

    iterated_solutions = iterated_solutions(:,1:iterations_to_converge);
    
    velocity_solutions(1,i) = iterated_solutions(1,length(iterated_solutions));
    velocity_solutions(2,i) = iterated_solutions(2,length(iterated_solutions));
    
    total_iterations = total_iterations + iterations_to_converge;
    
    clear iterated_solutions
    
    
    difference_between_solutions = abs( velocity_solutions(1,i) - velocity_solutions(2,i) );
    
    margin_of_criticality = 1; %generous for now
    if difference_between_solutions < margin_of_criticality
        %then it has gone critical in this place
        critical = true;
        critical_pressure_ratio_location = positions(i);
    end
    
    
end






%calculate the other variables


velocity_solutions_final = zeros(1,length(velocity_solutions));


for i = 1:length(positions)
    
    
   x = positions(i);
   A = (4/NozzleLength^2)*(AreaAtEnds-AreaAtThroat)*x^2 - (4/NozzleLength)*(AreaAtEnds-AreaAtThroat)*x + AreaAtEnds;
    
   
   if critical == true
       if x < critical_pressure_ratio_location
           velocity_solutions_final(i) = velocity_solutions(2,i);
       else
           velocity_solutions_final(i) = velocity_solutions(1,i);
       end
   else
       velocity_solutions_final(i) = velocity_solutions(2,i);
   end
   
   pressure_solutions(i) = p0^(1-g)*((T0*m*R)/(velocity_solutions_final(i)*A))^g;
   temperature_solutions(i) = velocity_solutions_final(i)*pressure_solutions(i)*A/(m*R);
   density_solutions(i) = pressure_solutions(i)/(R*temperature_solutions(i));
   mach_number_solutions(i) = velocity_solutions_final(i)/sqrt(g*R*temperature_solutions(i));
   speed_of_sound_solutions(i) = sqrt( g*pressure_solutions(i)/density_solutions(i) );
   area_ratio_solutions(i) = (1/mach_number_solutions(i)) * ( ((g+1)/2)/(1 + ((g-1)/2)*mach_number_solutions(i)^2) )^((g+1)/(2*(1-g)));
   
   
end






%0D capacity plot
figure(1)
hold on
title('Nozzle capacity vs pressure ratio')
ylabel('Capacity')
xlabel('Pressure ratio')
xlim([0,1])
set(gca,'FontSize',14)
xline(PressureRatioCritical,'--r');
yline(CapacityChoked,'--b');
plot(PressureRatiosArray,Capacities,'g','LineWidth',2);

%1D area vesus position
figure(2)
hold on
title('Nozzle area vs position')
ylabel('Area, m^2')
xlabel('Position')
set(gca,'FontSize',14)
plot(positions,Areas,'g','LineWidth',2);

%plot of position along nozzle vs velocity
figure(3)
hold on
title('Flow velocity vs position')
ylabel('Velocity, m/s')
xlabel('Position')
set(gca,'FontSize',14)
if critical == true
    xline(critical_pressure_ratio_location,'--r');
end
plot(positions,speed_of_sound_solutions,'--b');
plot(positions,velocity_solutions_final,'k','LineWidth',2);

%plot of position along nozzle vs mach number
figure(4)
hold on
title('Mach number vs position')
ylabel('Mach number')
xlabel('Position')
set(gca,'FontSize',14)
if critical == true
    xline(critical_pressure_ratio_location,'--r');
end
yline(1,'--b');
plot(positions,mach_number_solutions,'k','LineWidth',2);

%plot of position along nozzle vs pressure
figure(5)
hold on
title('Pressure ratio vs position')
ylabel('Pressure ratio')
xlabel('Position')
set(gca,'FontSize',14)
if critical == true
    xline(critical_pressure_ratio_location,'--r');
end
yline(PressureRatioCritical,'--b');
plot(positions,pressure_solutions./p0,'k','LineWidth',2);

%plot of position along nozzle vs pressure
figure(6)
hold on
title('Effective area ratio vs position')
ylabel('Ratio')
xlabel('Position')
set(gca,'FontSize',14)
%if critical == true
%    xline(critical_pressure_ratio_location,'--r');
%end
%yline(PressureRatioCritical,'--b');
plot(positions,area_ratio_solutions,'k','LineWidth',2);









