function [streamlineVelocities,streamlinePressures,streamlineTemperatures,streamlineDensities,streamlineMachNumbers,streamlineSpeedsOfSound,streamlineAreaRatios] = streamlineSolver_ver01(streamlineAreas,streamlineMassFlowRate)



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



velocity_solutions = zeros(2,length(streamlineAreas));
streamlinePressures = zeros(1,length(streamlineAreas));
streamlineDensities = zeros(1,length(streamlineAreas));
streamlineTemperatures = zeros(1,length(streamlineAreas));
streamlineMachNumbers = zeros(1,length(streamlineAreas));
streamlineSpeedsOfSound = zeros(1,length(streamlineAreas));
streamlineAreaRatios = zeros(1,length(streamlineAreas));



%track whether the streamline has gone critical
critical = false;
%maximum speed of the nozzle
starting_guess_high = sqrt(  2*p0/r0*(g/(g-1))  ); 
%minimum speed of the nozzle
starting_guess_low = 0;
%track the total number of iterations
total_iterations = 0;



for i = 1:length(streamlineAreas)
    
    
    %start at the starting guesses, but after the first answers, use each
    %velocity to guess the next one downstream
    iterated_solutions = zeros(2,1e6);
    if i == 1
        iterated_solutions(1,1) = starting_guess_high;
        iterated_solutions(2,1) = starting_guess_low;
    else
        iterated_solutions(1,1) = velocity_solutions(1,i-1);
        iterated_solutions(2,1) = velocity_solutions(2,i-1);
    end
    
    
    A = streamlineAreas(i);
    residual = starting_guess_high;
    
    n = 2;
    while residual > 1e-6
        
        %this equation is real between m/(r0*A) and infinity
        iterated_solutions(1,n) = sqrt(  2*(p0/r0)*(g/(g-1))*( 1 - (streamlineMassFlowRate/(r0*A*iterated_solutions(1,n-1)))^(g-1) )  );
        
        %this equation is real between 0 and the maximum nozzle speed
        
        iterated_solutions(2,n) = (streamlineMassFlowRate/(r0*A))*(  1 - 0.5*(r0/p0)*((g-1)/g)*iterated_solutions(2,n-1)^2  )^(1/(1-g)) ;
        residual = abs( iterated_solutions(2,n) - iterated_solutions(2,n-1) );
        
        n = n + 1;
    end
    
    iterations_to_converge = n-1;
    total_iterations = total_iterations + iterations_to_converge;
    
    iterated_solutions = iterated_solutions(:,1:iterations_to_converge);
    velocity_solutions(1,i) = iterated_solutions(1,length(iterated_solutions));
    velocity_solutions(2,i) = iterated_solutions(2,length(iterated_solutions));
    clear iterated_solutions
    
    
    difference_between_solutions = abs( velocity_solutions(1,i) - velocity_solutions(2,i) );
   
    margin_of_criticality = 1; %generous for now
    if difference_between_solutions < margin_of_criticality
        %then it has gone critical in this place
        critical = true;
        critical_pressure_ratio_location = i;
    end
    
    
end



%calculate the other variables


streamlineVelocities = zeros(1,length(velocity_solutions));


for i = 1:length(streamlineAreas)
    
    
   A = streamlineAreas(i);
   
   
   if critical == true
       if i < critical_pressure_ratio_location
           streamlineVelocities(i) = velocity_solutions(2,i);
       else
           streamlineVelocities(i) = velocity_solutions(1,i);
       end
   else
       streamlineVelocities(i) = velocity_solutions(2,i);
   end
   
   
   streamlinePressures(i) = p0^(1-g)*((T0*streamlineMassFlowRate*R)/(streamlineVelocities(i)*A))^g;
   streamlineTemperatures(i) = streamlineVelocities(i)*streamlinePressures(i)*A/(streamlineMassFlowRate*R);
   streamlineDensities(i) = streamlinePressures(i)/(R*streamlineTemperatures(i));
   streamlineMachNumbers(i) = streamlineVelocities(i)/sqrt(g*R*streamlineTemperatures(i));
   streamlineSpeedsOfSound(i) = sqrt( g*streamlinePressures(i)/streamlineDensities(i) );
   streamlineAreaRatios(i) = (1/streamlineMachNumbers(i)) * ( ((g+1)/2)/(1 + ((g-1)/2)*streamlineMachNumbers(i)^2) )^((g+1)/(2*(1-g)));
   
   
end



end