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
AreaAtThroat = 0.1;
AreaAtEnds = 0.2;









%Capacity vs pressure ratio for a simple nozzle


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





%1D area versus position





positions = 0:0.001:NozzleLength;

Areas = (4/NozzleLength^2)*(AreaAtEnds-AreaAtThroat).*positions.^2 - (4/NozzleLength)*(AreaAtEnds-AreaAtThroat).*positions + AreaAtEnds;
    








%%%%The 1D CFD

%throat pressure ratio
r = PressureRatioCritical;

%mass flow rate
m = (p0/sqrt(T0))*S*AreaAtThroat*sqrt(r^(2/g)*(1-r^((g-1)/g)));

velocity_solutions = zeros(2,length(positions));
pressure_solutions = zeros(2,length(positions));
density_solutions = zeros(2,length(positions));
temperature_solutions = zeros(2,length(positions));
mach_number_solutions = zeros(2,length(positions));
speed_of_sound_solutions = zeros(2,length(positions));


%initial range
possible_velocities = 0:1:1000; %solve numerically for the correct one

critical_pressure_ratio_location = 0;


for i = 1:length(positions)
    
    
    x = positions(i);
    A = (4/NozzleLength^2)*(AreaAtEnds-AreaAtThroat)*x^2 - (4/NozzleLength)*(AreaAtEnds-AreaAtThroat)*x + AreaAtEnds;
    
    
    this_solution_velocities = zeros(1,2);
    
    for j = 2:length(possible_velocities)
        
        v_last  = possible_velocities(j-1);
        v_this = possible_velocities(j);
       
        sign_last = sign(0.5*(g-1)/g*v_last^2 + R*T0*(m/(r0*A*v_last))^(g-1) - p0/r0);
        sign_this = sign(0.5*(g-1)/g*v_this^2 + R*T0*(m/(r0*A*v_this))^(g-1) - p0/r0);
        
        if sign_last ~= sign_this
            if this_solution_velocities(1) == 0
                this_solution_velocities(1) = 0.5*(v_last+v_this);
            else
                this_solution_velocities(2) = 0.5*(v_last+v_this);
            end
        end
        
    end
    
    
   velocity_solutions(1,i) = this_solution_velocities(1);
   velocity_solutions(2,i) = this_solution_velocities(2);
   
   
   %refine the range of possible velocities by assuming no sudden changes
   if this_solution_velocities(1) ~= 0 && this_solution_velocities(1) ~= 0
       bottom_of_range = this_solution_velocities(1)*0.9;
       top_of_range = this_solution_velocities(2)*1.1;
       grain_size = (top_of_range-bottom_of_range)/1000;
       possible_velocities = bottom_of_range:grain_size:top_of_range;
   else
       %then this is the critical pressure ratio
       critical_pressure_ratio_location = positions(i);
   end
    
   
end



%calculate the other variables


for i = 1:length(positions)
    
    
    %allow for the discontinuity at M=1
    if velocity_solutions(1,i) == 0 || velocity_solutions(2,i) == 0
        velocity_solutions(1,i) = 0.5*( velocity_solutions(1,i-1) + velocity_solutions(1,i+1) );
        velocity_solutions(2,i) = 0.5*( velocity_solutions(2,i-1) + velocity_solutions(2,i+1) );
    end

    
    x = positions(i);
    A = (4/NozzleLength^2)*(AreaAtEnds-AreaAtThroat)*x^2 - (4/NozzleLength)*(AreaAtEnds-AreaAtThroat)*x + AreaAtEnds;
    
    
   pressure_solutions(1,i) = p0^(1-g)*((T0*m*R)/(velocity_solutions(1,i)*A))^g;
   pressure_solutions(2,i) = p0^(1-g)*((T0*m*R)/(velocity_solutions(2,i)*A))^g;
   temperature_solutions(1,i) = velocity_solutions(1,i)*pressure_solutions(1,i)*A/(m*R);
   temperature_solutions(2,i) = velocity_solutions(2,i)*pressure_solutions(2,i)*A/(m*R);
   density_solutions(1,i) = pressure_solutions(1,i)/(R*temperature_solutions(1,i));
   density_solutions(2,i) = pressure_solutions(2,i)/(R*temperature_solutions(2,i));
   mach_number_solutions(1,i) = velocity_solutions(1,i)/sqrt(g*R*temperature_solutions(1,i));
   mach_number_solutions(2,i) = velocity_solutions(2,i)/sqrt(g*R*temperature_solutions(2,i));
   speed_of_sound_solutions(1,i) = sqrt( g*pressure_solutions(1,i)/density_solutions(1,i) );
   speed_of_sound_solutions(2,i) = sqrt( g*pressure_solutions(2,i)/density_solutions(2,i) );
   
   
end


%x = 0.5;
%A = (4/NozzleLength^2)*(AreaAtEnds-AreaAtThroat)*x^2 - (4/NozzleLength)*(AreaAtEnds-AreaAtThroat)*x + AreaAtEnds;

%test_velocities = 10:0.05:700; %solve numerically for the correct one
%test_output = 0.5*(g-1)/g.*test_velocities.^2 + R*T0*(m./(r0*A.*test_velocities)).^(g-1) - p0/r0;


%test_c = sqrt(2*(g/(g-1))*p0/r0);



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

%plot of position along nozzle vs velocity
figure(3)
hold on
if critical_pressure_ratio_location ~= 0
    xline(critical_pressure_ratio_location,'--r');
end
plot(positions,speed_of_sound_solutions,'--b');
plot(positions,velocity_solutions,'k','LineWidth',2);

%plot of position along nozzle vs mach number
figure(4)
hold on
if critical_pressure_ratio_location ~= 0
    xline(critical_pressure_ratio_location,'--r');
end
yline(1,'--b');
plot(positions,mach_number_solutions,'k','LineWidth',2);

%plot of position along nozzle vs pressure
figure(5)
hold on
yline(PressureRatioCritical,'--b');
plot(positions,pressure_solutions./p0,'k','LineWidth',2);

%testing both solutions
%figure(6)
%hold on
%yline(0,'--b');
%plot(test_velocities,left,'r','LineWidth',2);
%plot(test_velocities,right,'b','LineWidth',2);
%plot(test_velocities,test_output,'k','LineWidth',2);







