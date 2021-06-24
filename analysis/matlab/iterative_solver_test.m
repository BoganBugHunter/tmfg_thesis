



A = 0.1;

solutions = zeros(2,1e6);

%maximum speed of the nozzle
starting_guess_high = sqrt(  2*p0/r0*(g/(g-1))  ); 

%minimum speed of the nozzle
starting_guess_low = 0;

solutions(1,1) = starting_guess_high;
solutions(2,1) = starting_guess_low;

residual = starting_guess_high;

n = 2;

while residual > 1e-6
    
    %this equation is real between m/(r0*A) and infinity
    solutions(1,n) = sqrt(  2*(p0/r0)*(g/(g-1))*( 1 - (m/(r0*A*solutions(1,n-1)))^(g-1) )  );
    
    %this equation is real between 0 and the maximum nozzle speed
    solutions(2,n) = (m/(r0*A))*(  1 - 0.5*(r0/p0)*((g-1)/g)*solutions(2,n-1)^2  )^(1/(1-g)) ;
    
    residual = abs( solutions(2,n) - solutions(2,n-1) );
    
    solution_high = solutions(1,n);
    solution_low = solutions(2,n);
    
    n = n + 1;

end

iterations_to_converge = n - 1;

solutions = solutions(:,1:iterations_to_converge);





figure(1)
hold on
plot(solutions(1,:),'r','LineWidth',2);
plot(solutions(2,:),'b','LineWidth',2);

clear solutions


%minimum value before v-high becomes complex
%v_critical_high = m/(r0*A);

%maximum value before v-low becomes complex
%v_critical_low = sqrt(2*(p0/r0)*(g/(g-1)));


%test_input_1 = v_critical_high:1:8000;

%test1 = zeros(2,length(test_input_1));

%for n = 1:length(test_input_1)

%test1(1,n) = sqrt(  2*(p0/r0)*(g/(g-1))*( 1 - (m/(r0*A*test_input_1(n)))^(g-1) )  );

%end

%test_input_2 = 0:1:v_critical_low;

%test2 = zeros(2,length(test_input_2));

%for n = 1:length(test_input_2)

%test2(2,n) = (m/(r0*A))*(  1 - 0.5*(r0/p0)*((g-1)/g)*test_input_2(n)^2  )^(1/(1-g)) ;

%end


%figure(2)
%hold on
%plot(test_input_1,test1(1,:), 'r-')
%plot(test_input_2,test2(2,:), 'b-')
%plot(1:1:800,1:1:800, 'k-')








