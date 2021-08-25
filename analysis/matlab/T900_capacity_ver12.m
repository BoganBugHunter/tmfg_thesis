%T900_capacity
%plots the q-curves for the four T900 smooth vane midsections

%ver01 KSZ03
%ver02 gets rid of individual plots; all comparative
%ver03 calculates capacities and PRs before plotting, not during
%ver04 includes pseudo-3D approximation
%ver05 uses curve interpolation on capacity data and tidies up the plot commands
%ver06 has changed plots for R2014b and has the novel data
%ver07 has the 3D mean capacity scatter plots from digitised RR data
%ver08 has the realistic TE SS cutback study
%ver09 is lighter - unnecessary plots are gone
%ver10 has the TE blowing rate study
%ver11 has the energy loss analysis for variable cutbacks
%ver12 is the first version modified by the 2021 work

close all;
clear variables;
load workspaces/T900_rolls_ver17;

%T900_ver06b does not have the outlier in KSZ03
%T900_rolls_ver02 includes the banana-2D throat widths
%T900_rolls_ver03 includes the properly-2D throat widths
%T900_rolls_ver04 includes the 3D throat areas
%T900_rolls_ver05 fixes the extra row at PR = 1.6477
%T900_rolls_ver06 changes "T900_data_3D" to "T900_capacities_3D"
%T900_rolls_ver07 has the 3D data in seperate 2-column vectors
%T900_rolls_ver08 is corrected to include the one missing row of the 2D CFD data (now row 20) and has the novel TE with cutback 0
%T900_rolls_ver09 has novel TE with cutback 1 and includes the digitised RR design capacities
%T900_rolls_ver09b is playing around with the capacity data to see if
%moving them up or down by one PR step looks more correct and fixing the
%outlier in KSZ03 at 20. Conclusion: there was something a bit wrong with
%all of them. They have been re-entered from the raw data files and seem
%fine now
%T900_rolls_ver10 has adopted the changes from ver09b and also has cutback
%2
%T900_rolls_ver10b has the rerun of cutback 2 (T013) and cutback 3
%T900_rolls_ver10c has the remesh of cutback 4 (T016) in column 6 and the
%original cutback 4 (T015) in column 7 and the remesh of cutback 3 in
%column 8 and the remesh of cutback 2 in column 9 and the remesh of cutback
%1 in column 10
%T900_rolls_ver11 has the realistic TE SS cutback study
%T900_rolls_ver12 has the outlet loss analysis for the realistic TEs
%T900_rolls_ver13 has the unsteady data added to T900_TE_data_2D
%T900_rolls_ver14 has the TE blowing rate data
%T900_rolls_ver15 has the outlet static pressures at different cutbacks to
%do energy loss analysis
%T900_rolls_ver16 is the first workspace loaded by the 2021 work
%T900_rolls_ver17 includes the turning angles for the suction-side cutbacks
%and the total pressures on the domain outlet



T900_capacities_2D = T900_data_2D;
T900_capacities_2D(:,1) = 1./((T900_data_2D(:,1) + 4285900)./(39690 + 4285900));
T900_capacities_2D(:,2:end) = T900_data_2D(:,2:end)*sqrt(300)/(39690 + 4285900);

T900_novel_capacities_2D = T900_novel_data_2D;
T900_novel_capacities_2D(:,1) = 1./((T900_novel_data_2D(:,1) + 4285900)./(39690 + 4285900));
T900_novel_capacities_2D(:,2:end) = T900_novel_data_2D(:,2:end)*sqrt(300)/(39690 + 4285900);

%T900_TE_data_2D contains unsteady data from column 13 onwards
T900_TE_capacities_2D = T900_TE_data_2D;
T900_TE_capacities_2D(:,1) = 1./((T900_TE_data_2D(:,1) + 4285900)./(39690 + 4285900));
T900_TE_capacities_2D(:,2:end) = T900_TE_data_2D(:,2:end)*sqrt(300)/(39690 + 4285900);



%%%%%%%%%%%%%%%%%%%%%%MEAN CAPACITIES 2D%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Capacity/Area trends

T900_capacities_2D_design = T900_capacities_2D(22,:) + (0.05086/0.079022)*(T900_capacities_2D(22,:) - T900_capacities_2D(22,:));

T900_capacity_Mskew_mean = mean(T900_capacities_2D_design(2:4));
T900_capacity_EP1_mean = mean(T900_capacities_2D_design(5:7));

throat_width_Mskew_mean = mean(throat_widths(1:3));
throat_width_EP1_mean = mean(throat_widths(4:6));

%percentage changes
T900_capacity_change_mean = 100*((T900_capacity_EP1_mean - T900_capacity_Mskew_mean)/T900_capacity_Mskew_mean);
throat_width_change_mean = 100*((throat_width_EP1_mean - throat_width_Mskew_mean)/throat_width_Mskew_mean);

%throat number: capacity percentage change over width percentage change
throat_number_2D = T900_capacity_change_mean/throat_width_change_mean;



%%%%%%%%%%%%%%%%%%%%%%MEAN CAPACITIES 3D%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Capacity/Area trends

T900_capacity_Mskew_mean_3D = mean(T900_capacities_3D_design(2:4));
T900_capacity_EP1_mean_3D = mean(T900_capacities_3D_design(5:7));

throat_width_Mskew_mean_3D = mean(throat_areas_digitised(1:3));
throat_width_EP1_mean_3D = mean(throat_areas_digitised(4:6));

%percentage changes
T900_capacity_change_mean_3D = 100*((T900_capacity_EP1_mean_3D - T900_capacity_Mskew_mean_3D)/T900_capacity_Mskew_mean_3D);
throat_width_change_mean_3D = 100*((throat_width_EP1_mean_3D - throat_width_Mskew_mean_3D)/throat_width_Mskew_mean_3D);

%throat number: capacity percentage change over width percentage change
throat_number_3D = T900_capacity_change_mean_3D/throat_width_change_mean_3D;



%%%%%%%%%%%%OVERSAMPLING 2D%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


n = size(T900_capacities_2D);
m = n(2);
n = n(1);

range = T900_capacities_2D(n,1) - T900_capacities_2D(1,1);
  
resolution = 10000;
    
interp_vector = T900_capacities_2D(1,1):range/resolution:T900_capacities_2D(n,1);
interp_vector = interp_vector';

p = size(interp_vector);
p = p(1);

T900_capacities_2D_oversampled = zeros(p,m);

for i2 = 1:m

T900_capacities_2D_oversampled(:,i2) = interp1(T900_capacities_2D(:,1),T900_capacities_2D(:,i2),interp_vector, 'pchip');

end  


n = size(T900_novel_capacities_2D);
m = n(2);
n = n(1);

range = T900_novel_capacities_2D(n,1) - T900_novel_capacities_2D(1,1);
  
resolution = 10000;
    
interp_vector = T900_novel_capacities_2D(1,1):range/resolution:T900_novel_capacities_2D(n,1);
interp_vector = interp_vector';

p = size(interp_vector);
p = p(1);

T900_novel_capacities_2D_oversampled = zeros(p,m);

for i2 = 1:m

T900_novel_capacities_2D_oversampled(:,i2) = interp1(T900_novel_capacities_2D(:,1),T900_novel_capacities_2D(:,i2),interp_vector, 'pchip');

end  


n = size(T900_TE_capacities_2D);
m = n(2);
n = n(1);


range = T900_TE_capacities_2D(n,1) - T900_TE_capacities_2D(1,1);
  
resolution = 10000;
    
interp_vector = T900_TE_capacities_2D(1,1):range/resolution:T900_TE_capacities_2D(n,1);
interp_vector = interp_vector';

p = size(interp_vector);
p = p(1);

T900_TE_capacities_2D_oversampled = zeros(p,m);

for i2 = 1:m

T900_TE_capacities_2D_oversampled(:,i2) = interp1(T900_TE_capacities_2D(:,1),T900_TE_capacities_2D(:,i2),interp_vector, 'pchip');

end 



%%%%%%%%%%%OVERSAMPLING 3D%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


resolution = 1000;


n = size(T900_capacity_3D_1);
m = n(2);
n = n(1);

range = T900_capacity_3D_1(n,1) - T900_capacity_3D_1(1,1);
    
interp_vector = T900_capacity_3D_1(1,1):range/resolution:T900_capacity_3D_1(n,1);
interp_vector = interp_vector';

p = size(interp_vector);
p = p(1);

T900_capacity_3D_1_oversampled = zeros(p,m);

for i2 = 1:m

T900_capacity_3D_1_oversampled(:,i2) = interp1(T900_capacity_3D_1(:,1),T900_capacity_3D_1(:,i2),interp_vector, 'pchip');

end  


n = size(T900_capacity_3D_2);
m = n(2);
n = n(1);

range = T900_capacity_3D_2(n,1) - T900_capacity_3D_2(1,1);
    
interp_vector = T900_capacity_3D_2(1,1):range/resolution:T900_capacity_3D_2(n,1);
interp_vector = interp_vector';

p = size(interp_vector);
p = p(1);

T900_capacity_3D_2_oversampled = zeros(p,m);

for i2 = 1:m

T900_capacity_3D_2_oversampled(:,i2) = interp1(T900_capacity_3D_2(:,1),T900_capacity_3D_2(:,i2),interp_vector, 'pchip');

end  


n = size(T900_capacity_3D_3);
m = n(2);
n = n(1);

range = T900_capacity_3D_3(n,1) - T900_capacity_3D_3(1,1);
    
interp_vector = T900_capacity_3D_3(1,1):range/resolution:T900_capacity_3D_3(n,1);
interp_vector = interp_vector';

p = size(interp_vector);
p = p(1);

T900_capacity_3D_3_oversampled = zeros(p,m);

for i2 = 1:m

T900_capacity_3D_3_oversampled(:,i2) = interp1(T900_capacity_3D_3(:,1),T900_capacity_3D_3(:,i2),interp_vector, 'pchip');

end  


n = size(T900_capacity_3D_4);
m = n(2);
n = n(1);

range = T900_capacity_3D_4(n,1) - T900_capacity_3D_4(1,1);
    
interp_vector = T900_capacity_3D_4(1,1):range/resolution:T900_capacity_3D_4(n,1);
interp_vector = interp_vector';

p = size(interp_vector);
p = p(1);

T900_capacity_3D_4_oversampled = zeros(p,m);

for i2 = 1:m

T900_capacity_3D_4_oversampled(:,i2) = interp1(T900_capacity_3D_4(:,1),T900_capacity_3D_4(:,i2),interp_vector, 'pchip');

end  


n = size(T900_capacity_3D_5);
m = n(2);
n = n(1);

range = T900_capacity_3D_5(n,1) - T900_capacity_3D_5(1,1);
    
interp_vector = T900_capacity_3D_5(1,1):range/resolution:T900_capacity_3D_5(n,1);
interp_vector = interp_vector';

p = size(interp_vector);
p = p(1);

T900_capacity_3D_5_oversampled = zeros(p,m);

for i2 = 1:m

T900_capacity_3D_5_oversampled(:,i2) = interp1(T900_capacity_3D_5(:,1),T900_capacity_3D_5(:,i2),interp_vector, 'pchip');

end  


n = size(T900_capacity_3D_6);
m = n(2);
n = n(1);

range = T900_capacity_3D_6(n,1) - T900_capacity_3D_6(1,1);
    
interp_vector = T900_capacity_3D_6(1,1):range/resolution:T900_capacity_3D_6(n,1);
interp_vector = interp_vector';

p = size(interp_vector);
p = p(1);

T900_capacity_3D_6_oversampled = zeros(p,m);

for i2 = 1:m

T900_capacity_3D_6_oversampled(:,i2) = interp1(T900_capacity_3D_6(:,1),T900_capacity_3D_6(:,i2),interp_vector, 'pchip');

end 



%%%%%%%%%%%PSEUDO-3D_APPROXIMATION%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%number of points to look at on Q-curve. range is "range_on_x_axis" above and below the central point
range_on_x_axis = 2000;

T900_capacities_P3D = T900_capacities_2D_oversampled;

n = size(T900_capacities_2D_oversampled);


for i1 = (range_on_x_axis+1):(n-range_on_x_axis-1)
    
    sample = T900_capacities_2D_oversampled(i1-range_on_x_axis:i1+range_on_x_axis,2:7);
    
    T900_capacities_P3D(i1,2:7) = sum(sample,1)/(2*range_on_x_axis + 1);

end



%%%%%%%%%%%%VARIABLE CUTBACKS - TOTAL PRESSURE INTEGRATORS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


distances_10 = outlet_totalpressure_10(:,1);
pressures_10 = outlet_totalpressure_10(:,2);

array_size = size(distances_10);
array_size = array_size(1);

integrand_array = zeros(array_size,1);

for n = 1:(array_size-1)

integrand_array(n) = 0.5*(distances_10(n+1) - distances_10(n))*(pressures_10(n+1) + pressures_10(n));

end

total_pressure_integrated_10 = sum(integrand_array);


distances_9 = outlet_totalpressure_9(:,1);
pressures_9 = outlet_totalpressure_9(:,2);

array_size = size(distances_9);
array_size = array_size(1);

integrand_array = zeros(array_size,1);

for n = 1:(array_size-1)

integrand_array(n) = 0.5*(distances_9(n+1) - distances_9(n))*(pressures_9(n+1) + pressures_9(n));

end

total_pressure_integrated_9 = sum(integrand_array);


distances_8 = outlet_totalpressure_8(:,1);
pressures_8 = outlet_totalpressure_8(:,2);

array_size = size(distances_8);
array_size = array_size(1);

integrand_array = zeros(array_size,1);

for n = 1:(array_size-1)

integrand_array(n) = 0.5*(distances_8(n+1) - distances_8(n))*(pressures_8(n+1) + pressures_8(n));

end

total_pressure_integrated_8 = sum(integrand_array);


distances_7 = outlet_totalpressure_7(:,1);
pressures_7 = outlet_totalpressure_7(:,2);

array_size = size(distances_7);
array_size = array_size(1);

integrand_array = zeros(array_size,1);

for n = 1:(array_size-1)

integrand_array(n) = 0.5*(distances_7(n+1) - distances_7(n))*(pressures_7(n+1) + pressures_7(n));

end

total_pressure_integrated_7 = sum(integrand_array);


distances_6 = outlet_totalpressure_6(:,1);
pressures_6 = outlet_totalpressure_6(:,2);

array_size = size(distances_6);
array_size = array_size(1);

integrand_array = zeros(array_size,1);

for n = 1:(array_size-1)

integrand_array(n) = 0.5*(distances_6(n+1) - distances_6(n))*(pressures_6(n+1) + pressures_6(n));

end

total_pressure_integrated_6 = sum(integrand_array);


distances_5 = outlet_totalpressure_5(:,1);
pressures_5 = outlet_totalpressure_5(:,2);

array_size = size(distances_5);
array_size = array_size(1);

integrand_array = zeros(array_size,1);

for n = 1:(array_size-1)

integrand_array(n) = 0.5*(distances_5(n+1) - distances_5(n))*(pressures_5(n+1) + pressures_5(n));

end

total_pressure_integrated_5 = sum(integrand_array);


distances_4 = outlet_totalpressure_4(:,1);
pressures_4 = outlet_totalpressure_4(:,2);

array_size = size(distances_4);
array_size = array_size(1);

integrand_array = zeros(array_size,1);

for n = 1:(array_size-1)

integrand_array(n) = 0.5*(distances_4(n+1) - distances_4(n))*(pressures_4(n+1) + pressures_4(n));

end

total_pressure_integrated_4 = sum(integrand_array);


distances_3 = outlet_totalpressure_3(:,1);
pressures_3 = outlet_totalpressure_3(:,2);

array_size = size(distances_3);
array_size = array_size(1);

integrand_array = zeros(array_size,1);

for n = 1:(array_size-1)

integrand_array(n) = 0.5*(distances_3(n+1) - distances_3(n))*(pressures_3(n+1) + pressures_3(n));

end

total_pressure_integrated_3 = sum(integrand_array);


distances_2 = outlet_totalpressure_2(:,1);
pressures_2 = outlet_totalpressure_2(:,2);

array_size = size(distances_2);
array_size = array_size(1);

integrand_array = zeros(array_size,1);

for n = 1:(array_size-1)

integrand_array(n) = 0.5*(distances_2(n+1) - distances_2(n))*(pressures_2(n+1) + pressures_2(n));

end

total_pressure_integrated_2 = sum(integrand_array);


distances_1 = outlet_totalpressure_1(:,1);
pressures_1 = outlet_totalpressure_1(:,2);

array_size = size(distances_1);
array_size = array_size(1);

integrand_array = zeros(array_size,1);

for n = 1:(array_size-1)

integrand_array(n) = 0.5*(distances_1(n+1) - distances_1(n))*(pressures_1(n+1) + pressures_1(n));

end

total_pressure_integrated_1 = sum(integrand_array);


distances_0 = outlet_totalpressure_0(:,1);
pressures_0 = outlet_totalpressure_0(:,2);

array_size = size(distances_0);
array_size = array_size(1);

integrand_array = zeros(array_size,1);

for n = 1:(array_size-1)

integrand_array(n) = 0.5*(distances_0(n+1) - distances_0(n))*(pressures_0(n+1) + pressures_0(n));

end

total_pressure_integrated_0 = sum(integrand_array);


distances_10_design = outlet_totalpressure_10_design(:,1);
pressures_10_design = outlet_totalpressure_10_design(:,2);

array_size = size(distances_10_design);
array_size = array_size(1);

integrand_array = zeros(array_size,1);

for n = 1:(array_size-1)

integrand_array(n) = 0.5*(distances_10_design(n+1) - distances_10_design(n))*(pressures_10_design(n+1) + pressures_10_design(n));

end

total_pressure_integrated_10_design = sum(integrand_array);


distances_9_design = outlet_totalpressure_9_design(:,1);
pressures_9_design = outlet_totalpressure_9_design(:,2);

array_size = size(distances_9_design);
array_size = array_size(1);

integrand_array = zeros(array_size,1);

for n = 1:(array_size-1)

integrand_array(n) = 0.5*(distances_9_design(n+1) - distances_9_design(n))*(pressures_9_design(n+1) + pressures_9_design(n));

end

total_pressure_integrated_9_design = sum(integrand_array);


distances_8_design = outlet_totalpressure_8_design(:,1);
pressures_8_design = outlet_totalpressure_8_design(:,2);

array_size = size(distances_8_design);
array_size = array_size(1);

integrand_array = zeros(array_size,1);

for n = 1:(array_size-1)

integrand_array(n) = 0.5*(distances_8_design(n+1) - distances_8_design(n))*(pressures_8_design(n+1) + pressures_8_design(n));

end

total_pressure_integrated_8_design = sum(integrand_array);


distances_7_design = outlet_totalpressure_7_design(:,1);
pressures_7_design = outlet_totalpressure_7_design(:,2);

array_size = size(distances_7_design);
array_size = array_size(1);

integrand_array = zeros(array_size,1);

for n = 1:(array_size-1)

integrand_array(n) = 0.5*(distances_7_design(n+1) - distances_7_design(n))*(pressures_7_design(n+1) + pressures_7_design(n));

end

total_pressure_integrated_7_design = sum(integrand_array);


distances_6_design = outlet_totalpressure_6_design(:,1);
pressures_6_design = outlet_totalpressure_6_design(:,2);

array_size = size(distances_6_design);
array_size = array_size(1);

integrand_array = zeros(array_size,1);

for n = 1:(array_size-1)

integrand_array(n) = 0.5*(distances_6_design(n+1) - distances_6_design(n))*(pressures_6_design(n+1) + pressures_6_design(n));

end

total_pressure_integrated_6_design = sum(integrand_array);


distances_5_design = outlet_totalpressure_5_design(:,1);
pressures_5_design = outlet_totalpressure_5_design(:,2);

array_size = size(distances_5_design);
array_size = array_size(1);

integrand_array = zeros(array_size,1);

for n = 1:(array_size-1)

integrand_array(n) = 0.5*(distances_5_design(n+1) - distances_5_design(n))*(pressures_5_design(n+1) + pressures_5_design(n));

end

total_pressure_integrated_5_design = sum(integrand_array);


distances_4_design = outlet_totalpressure_4_design(:,1);
pressures_4_design = outlet_totalpressure_4_design(:,2);

array_size = size(distances_4_design);
array_size = array_size(1);

integrand_array = zeros(array_size,1);

for n = 1:(array_size-1)

integrand_array(n) = 0.5*(distances_4_design(n+1) - distances_4_design(n))*(pressures_4_design(n+1) + pressures_4_design(n));

end

total_pressure_integrated_4_design = sum(integrand_array);


distances_3_design = outlet_totalpressure_3_design(:,1);
pressures_3_design = outlet_totalpressure_3_design(:,2);

array_size = size(distances_3_design);
array_size = array_size(1);

integrand_array = zeros(array_size,1);

for n = 1:(array_size-1)

integrand_array(n) = 0.5*(distances_3_design(n+1) - distances_3_design(n))*(pressures_3_design(n+1) + pressures_3_design(n));

end

total_pressure_integrated_3_design = sum(integrand_array);


distances_2_design = outlet_totalpressure_2_design(:,1);
pressures_2_design = outlet_totalpressure_2_design(:,2);

array_size = size(distances_2_design);
array_size = array_size(1);

integrand_array = zeros(array_size,1);

for n = 1:(array_size-1)

integrand_array(n) = 0.5*(distances_2_design(n+1) - distances_2_design(n))*(pressures_2_design(n+1) + pressures_2_design(n));

end

total_pressure_integrated_2_design = sum(integrand_array);


distances_1_design = outlet_totalpressure_1_design(:,1);
pressures_1_design = outlet_totalpressure_1_design(:,2);

array_size = size(distances_1_design);
array_size = array_size(1);

integrand_array = zeros(array_size,1);

for n = 1:(array_size-1)

integrand_array(n) = 0.5*(distances_1_design(n+1) - distances_1_design(n))*(pressures_1_design(n+1) + pressures_1_design(n));

end

total_pressure_integrated_1_design = sum(integrand_array);


distances_0_design = outlet_totalpressure_0_design(:,1);
pressures_0_design = outlet_totalpressure_0_design(:,2);

array_size = size(distances_0_design);
array_size = array_size(1);

integrand_array = zeros(array_size,1);

for n = 1:(array_size-1)

integrand_array(n) = 0.5*(distances_0_design(n+1) - distances_0_design(n))*(pressures_0_design(n+1) + pressures_0_design(n));

end

total_pressure_integrated_0_design = sum(integrand_array);




%%%%%%%%%%%%VARIABLE CUTBACKS - STATIC PRESSURE INTEGRATORS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


distances_10_design = outlet_staticpressure_10_design(:,1);
pressures_10_design = outlet_staticpressure_10_design(:,2);

array_size = size(distances_10_design);
array_size = array_size(1);

integrand_array = zeros(array_size,1);

for n = 1:(array_size-1)

integrand_array(n) = 0.5*(distances_10_design(n+1) - distances_10_design(n))*(pressures_10_design(n+1) + pressures_10_design(n));

end

static_pressure_integrated_10_design = sum(integrand_array);


distances_9_design = outlet_staticpressure_9_design(:,1);
pressures_9_design = outlet_staticpressure_9_design(:,2);

array_size = size(distances_9_design);
array_size = array_size(1);

integrand_array = zeros(array_size,1);

for n = 1:(array_size-1)

integrand_array(n) = 0.5*(distances_9_design(n+1) - distances_9_design(n))*(pressures_9_design(n+1) + pressures_9_design(n));

end

static_pressure_integrated_9_design = sum(integrand_array);


distances_8_design = outlet_staticpressure_8_design(:,1);
pressures_8_design = outlet_staticpressure_8_design(:,2);

array_size = size(distances_8_design);
array_size = array_size(1);

integrand_array = zeros(array_size,1);

for n = 1:(array_size-1)

integrand_array(n) = 0.5*(distances_8_design(n+1) - distances_8_design(n))*(pressures_8_design(n+1) + pressures_8_design(n));

end

static_pressure_integrated_8_design = sum(integrand_array);


distances_7_design = outlet_staticpressure_7_design(:,1);
pressures_7_design = outlet_staticpressure_7_design(:,2);

array_size = size(distances_7_design);
array_size = array_size(1);

integrand_array = zeros(array_size,1);

for n = 1:(array_size-1)

integrand_array(n) = 0.5*(distances_7_design(n+1) - distances_7_design(n))*(pressures_7_design(n+1) + pressures_7_design(n));

end

static_pressure_integrated_7_design = sum(integrand_array);


distances_6_design = outlet_staticpressure_6_design(:,1);
pressures_6_design = outlet_staticpressure_6_design(:,2);

array_size = size(distances_6_design);
array_size = array_size(1);

integrand_array = zeros(array_size,1);

for n = 1:(array_size-1)

integrand_array(n) = 0.5*(distances_6_design(n+1) - distances_6_design(n))*(pressures_6_design(n+1) + pressures_6_design(n));

end

static_pressure_integrated_6_design = sum(integrand_array);


distances_5_design = outlet_staticpressure_5_design(:,1);
pressures_5_design = outlet_staticpressure_5_design(:,2);

array_size = size(distances_5_design);
array_size = array_size(1);

integrand_array = zeros(array_size,1);

for n = 1:(array_size-1)

integrand_array(n) = 0.5*(distances_5_design(n+1) - distances_5_design(n))*(pressures_5_design(n+1) + pressures_5_design(n));

end

static_pressure_integrated_5_design = sum(integrand_array);


distances_4_design = outlet_staticpressure_4_design(:,1);
pressures_4_design = outlet_staticpressure_4_design(:,2);

array_size = size(distances_4_design);
array_size = array_size(1);

integrand_array = zeros(array_size,1);

for n = 1:(array_size-1)

integrand_array(n) = 0.5*(distances_4_design(n+1) - distances_4_design(n))*(pressures_4_design(n+1) + pressures_4_design(n));

end

static_pressure_integrated_4_design = sum(integrand_array);


distances_3_design = outlet_staticpressure_3_design(:,1);
pressures_3_design = outlet_staticpressure_3_design(:,2);

array_size = size(distances_3_design);
array_size = array_size(1);

integrand_array = zeros(array_size,1);

for n = 1:(array_size-1)

integrand_array(n) = 0.5*(distances_3_design(n+1) - distances_3_design(n))*(pressures_3_design(n+1) + pressures_3_design(n));

end

static_pressure_integrated_3_design = sum(integrand_array);


distances_2_design = outlet_staticpressure_2_design(:,1);
pressures_2_design = outlet_staticpressure_2_design(:,2);

array_size = size(distances_2_design);
array_size = array_size(1);

integrand_array = zeros(array_size,1);

for n = 1:(array_size-1)

integrand_array(n) = 0.5*(distances_2_design(n+1) - distances_2_design(n))*(pressures_2_design(n+1) + pressures_2_design(n));

end

static_pressure_integrated_2_design = sum(integrand_array);


distances_1_design = outlet_staticpressure_1_design(:,1);
pressures_1_design = outlet_staticpressure_1_design(:,2);

array_size = size(distances_1_design);
array_size = array_size(1);

integrand_array = zeros(array_size,1);

for n = 1:(array_size-1)

integrand_array(n) = 0.5*(distances_1_design(n+1) - distances_1_design(n))*(pressures_1_design(n+1) + pressures_1_design(n));

end

static_pressure_integrated_1_design = sum(integrand_array);


distances_0_design = outlet_staticpressure_0_design(:,1);
pressures_0_design = outlet_staticpressure_0_design(:,2);

array_size = size(distances_0_design);
array_size = array_size(1);

integrand_array = zeros(array_size,1);

for n = 1:(array_size-1)

integrand_array(n) = 0.5*(distances_0_design(n+1) - distances_0_design(n))*(pressures_0_design(n+1) + pressures_0_design(n));

end

static_pressure_integrated_0_design = sum(integrand_array);






%%%%%%%%%%%COMBINING THE LOSSES INTO ONE ARRAY%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


total_pressures_integrated(1) = total_pressure_integrated_0;
total_pressures_integrated(2) = total_pressure_integrated_1;
total_pressures_integrated(3) = total_pressure_integrated_2;
total_pressures_integrated(4) = total_pressure_integrated_3;
total_pressures_integrated(5) = total_pressure_integrated_4;
total_pressures_integrated(6) = total_pressure_integrated_5;
total_pressures_integrated(7) = total_pressure_integrated_6;
total_pressures_integrated(8) = total_pressure_integrated_7;
total_pressures_integrated(9) = total_pressure_integrated_8;
total_pressures_integrated(10) = total_pressure_integrated_9;
total_pressures_integrated(11) = total_pressure_integrated_10;

total_pressures_integrated_design(1) = total_pressure_integrated_0_design;
total_pressures_integrated_design(2) = total_pressure_integrated_1_design;
total_pressures_integrated_design(3) = total_pressure_integrated_2_design;
total_pressures_integrated_design(4) = total_pressure_integrated_3_design;
total_pressures_integrated_design(5) = total_pressure_integrated_4_design;
total_pressures_integrated_design(6) = total_pressure_integrated_5_design;
total_pressures_integrated_design(7) = total_pressure_integrated_6_design;
total_pressures_integrated_design(8) = total_pressure_integrated_7_design;
total_pressures_integrated_design(9) = total_pressure_integrated_8_design;
total_pressures_integrated_design(10) = total_pressure_integrated_9_design;
total_pressures_integrated_design(11) = total_pressure_integrated_10_design;

static_pressures_integrated_design(1) = static_pressure_integrated_0_design;
static_pressures_integrated_design(2) = static_pressure_integrated_1_design;
static_pressures_integrated_design(3) = static_pressure_integrated_2_design;
static_pressures_integrated_design(4) = static_pressure_integrated_3_design;
static_pressures_integrated_design(5) = static_pressure_integrated_4_design;
static_pressures_integrated_design(6) = static_pressure_integrated_5_design;
static_pressures_integrated_design(7) = static_pressure_integrated_6_design;
static_pressures_integrated_design(8) = static_pressure_integrated_7_design;
static_pressures_integrated_design(9) = static_pressure_integrated_8_design;
static_pressures_integrated_design(10) = static_pressure_integrated_9_design;
static_pressures_integrated_design(11) = static_pressure_integrated_10_design;



loss_array_size = size(total_pressures_integrated);
loss_array_size_design = size(total_pressures_integrated_design);

energy_loss_array_size_design = size(total_pressures_integrated_design);

loss_array = 100*(ones(loss_array_size) - total_pressures_integrated./(4325590*0.0588106)) - 100;
%4325590*0.0588106 is the upstream TP times the inlet width
loss_array_design = 100*(ones(loss_array_size_design) - total_pressures_integrated_design./(4325590*0.0588106)) - 100;


%energy_loss_array = (0.01*100)*(ones(loss_array_size) - ((ones(loss_array_size) - ((static_pressures_integrated)./(total_pressures_integrated)).^0.285714)./(ones(loss_array_size) - ((static_pressures_integrated)./(4325590)).^0.285714)));
%needs a fudge factor of 0.01 for some reason

energy_loss_array_design = (0.01*100)*(ones(energy_loss_array_size_design) - ((ones(energy_loss_array_size_design) - ((static_pressures_integrated_design)./(total_pressures_integrated_design)).^0.285714)./(ones(energy_loss_array_size_design) - ((static_pressures_integrated_design)./(4325590)).^0.285714)));
%needs a fudge factor of 0.01 for some reason




%%%%%%%%%%%%%FILTERING THE VANE STATIC PRESSURE HISTORY%%%%%%%%%%%%%%%%%%%%


%[B,A] = butter(15,0.09,'high');
%[z,p,k] = butter(15,0.09,'high');
%sos = zp2sos(z,p,k);
%fvtool(sos, 'Analysis', 'freq')
%vane_static_pressure_history_2(:,3) = filter(B,A,vane_static_pressure_history_2(:,2)+178000);
%tapping_static_pressure_history(:,3) = filter(B,A,tapping_static_pressure_history(:,2)+178000);





%%%%%%%%%%%%%%LOSS VS TE BLOWING STUDY - TOTAL PRESSURE INTEGRATORS%%%%%%%%


distances_blowing_0 = outlet_totalpressure_blowing_0(:,1);
pressures_blowing_0 = outlet_totalpressure_blowing_0(:,2);

array_size = size(distances_blowing_0);
array_size = array_size(1);

integrand_array = zeros(array_size,1);

for n = 1:(array_size-1)

integrand_array(n) = 0.5*(distances_blowing_0(n+1) - distances_blowing_0(n))*(pressures_blowing_0(n+1) + pressures_blowing_0(n));

end

total_pressure_integrated_blowing_0 = sum(integrand_array);


distances_blowing_0pt5 = outlet_totalpressure_blowing_0pt5(:,1);
pressures_blowing_0pt5 = outlet_totalpressure_blowing_0pt5(:,2);

array_size = size(distances_blowing_0pt5);
array_size = array_size(1);

integrand_array = zeros(array_size,1);

for n = 1:(array_size-1)

integrand_array(n) = 0.5*(distances_blowing_0pt5(n+1) - distances_blowing_0pt5(n))*(pressures_blowing_0pt5(n+1) + pressures_blowing_0pt5(n));

end

total_pressure_integrated_blowing_0pt5 = sum(integrand_array);


distances_blowing_1 = outlet_totalpressure_blowing_1(:,1);
pressures_blowing_1 = outlet_totalpressure_blowing_1(:,2);

array_size = size(distances_blowing_1);
array_size = array_size(1);

integrand_array = zeros(array_size,1);

for n = 1:(array_size-1)

integrand_array(n) = 0.5*(distances_blowing_1(n+1) - distances_blowing_1(n))*(pressures_blowing_1(n+1) + pressures_blowing_1(n));

end

total_pressure_integrated_blowing_1 = sum(integrand_array);


distances_blowing_1pt5 = outlet_totalpressure_blowing_1pt5(:,1);
pressures_blowing_1pt5 = outlet_totalpressure_blowing_1pt5(:,2);

array_size = size(distances_blowing_1pt5);
array_size = array_size(1);

integrand_array = zeros(array_size,1);

for n = 1:(array_size-1)

integrand_array(n) = 0.5*(distances_blowing_1pt5(n+1) - distances_blowing_1pt5(n))*(pressures_blowing_1pt5(n+1) + pressures_blowing_1pt5(n));

end

total_pressure_integrated_blowing_1pt5 = sum(integrand_array);


distances_blowing_2 = outlet_totalpressure_blowing_2(:,1);
pressures_blowing_2 = outlet_totalpressure_blowing_2(:,2);

array_size = size(distances_blowing_2);
array_size = array_size(1);

integrand_array = zeros(array_size,1);

for n = 1:(array_size-1)

integrand_array(n) = 0.5*(distances_blowing_2(n+1) - distances_blowing_2(n))*(pressures_blowing_2(n+1) + pressures_blowing_2(n));

end

total_pressure_integrated_blowing_2 = sum(integrand_array);


distances_blowing_2pt5 = outlet_totalpressure_blowing_2pt5(:,1);
pressures_blowing_2pt5 = outlet_totalpressure_blowing_2pt5(:,2);

array_size = size(distances_blowing_2pt5);
array_size = array_size(1);

integrand_array = zeros(array_size,1);

for n = 1:(array_size-1)

integrand_array(n) = 0.5*(distances_blowing_2pt5(n+1) - distances_blowing_2pt5(n))*(pressures_blowing_2pt5(n+1) + pressures_blowing_2pt5(n));

end

total_pressure_integrated_blowing_2pt5 = sum(integrand_array);


distances_blowing_3 = outlet_totalpressure_blowing_3(:,1);
pressures_blowing_3 = outlet_totalpressure_blowing_3(:,2);

array_size = size(distances_blowing_3);
array_size = array_size(1);

integrand_array = zeros(array_size,1);

for n = 1:(array_size-1)

integrand_array(n) = 0.5*(distances_blowing_3(n+1) - distances_blowing_3(n))*(pressures_blowing_3(n+1) + pressures_blowing_3(n));

end

total_pressure_integrated_blowing_3 = sum(integrand_array);


distances_blowing_3pt5 = outlet_totalpressure_blowing_3pt5(:,1);
pressures_blowing_3pt5 = outlet_totalpressure_blowing_3pt5(:,2);

array_size = size(distances_blowing_3pt5);
array_size = array_size(1);

integrand_array = zeros(array_size,1);

for n = 1:(array_size-1)

integrand_array(n) = 0.5*(distances_blowing_3pt5(n+1) - distances_blowing_3pt5(n))*(pressures_blowing_3pt5(n+1) + pressures_blowing_3pt5(n));

end

total_pressure_integrated_blowing_3pt5 = sum(integrand_array);


distances_blowing_4 = outlet_totalpressure_blowing_4(:,1);
pressures_blowing_4 = outlet_totalpressure_blowing_4(:,2);

array_size = size(distances_blowing_4);
array_size = array_size(1);

integrand_array = zeros(array_size,1);

for n = 1:(array_size-1)

integrand_array(n) = 0.5*(distances_blowing_4(n+1) - distances_blowing_4(n))*(pressures_blowing_4(n+1) + pressures_blowing_4(n));

end

total_pressure_integrated_blowing_4 = sum(integrand_array);



%%%%%%%%%%%%%%LOSS VS TE BLOWING STUDY - STATIC PRESSURE INTEGRATORS%%%%%%%%


distances_blowing_0 = outlet_staticpressure_blowing_0(:,1);
pressures_blowing_0 = outlet_staticpressure_blowing_0(:,2);

array_size = size(distances_blowing_0);
array_size = array_size(1);

integrand_array = zeros(array_size,1);

for n = 1:(array_size-1)

integrand_array(n) = 0.5*(distances_blowing_0(n+1) - distances_blowing_0(n))*(pressures_blowing_0(n+1) + pressures_blowing_0(n));

end

static_pressure_integrated_blowing_0 = sum(integrand_array);


distances_blowing_0pt5 = outlet_staticpressure_blowing_0pt5(:,1);
pressures_blowing_0pt5 = outlet_staticpressure_blowing_0pt5(:,2);

array_size = size(distances_blowing_0pt5);
array_size = array_size(1);

integrand_array = zeros(array_size,1);

for n = 1:(array_size-1)

integrand_array(n) = 0.5*(distances_blowing_0pt5(n+1) - distances_blowing_0pt5(n))*(pressures_blowing_0pt5(n+1) + pressures_blowing_0pt5(n));

end

static_pressure_integrated_blowing_0pt5 = sum(integrand_array);


distances_blowing_1 = outlet_staticpressure_blowing_1(:,1);
pressures_blowing_1 = outlet_staticpressure_blowing_1(:,2);

array_size = size(distances_blowing_1);
array_size = array_size(1);

integrand_array = zeros(array_size,1);

for n = 1:(array_size-1)

integrand_array(n) = 0.5*(distances_blowing_1(n+1) - distances_blowing_1(n))*(pressures_blowing_1(n+1) + pressures_blowing_1(n));

end

static_pressure_integrated_blowing_1 = sum(integrand_array);


distances_blowing_1pt5 = outlet_staticpressure_blowing_1pt5(:,1);
pressures_blowing_1pt5 = outlet_staticpressure_blowing_1pt5(:,2);

array_size = size(distances_blowing_1pt5);
array_size = array_size(1);

integrand_array = zeros(array_size,1);

for n = 1:(array_size-1)

integrand_array(n) = 0.5*(distances_blowing_1pt5(n+1) - distances_blowing_1pt5(n))*(pressures_blowing_1pt5(n+1) + pressures_blowing_1pt5(n));

end

static_pressure_integrated_blowing_1pt5 = sum(integrand_array);


distances_blowing_2 = outlet_staticpressure_blowing_2(:,1);
pressures_blowing_2 = outlet_staticpressure_blowing_2(:,2);

array_size = size(distances_blowing_2);
array_size = array_size(1);

integrand_array = zeros(array_size,1);

for n = 1:(array_size-1)

integrand_array(n) = 0.5*(distances_blowing_2(n+1) - distances_blowing_2(n))*(pressures_blowing_2(n+1) + pressures_blowing_2(n));

end

static_pressure_integrated_blowing_2 = sum(integrand_array);


distances_blowing_2pt5 = outlet_staticpressure_blowing_2pt5(:,1);
pressures_blowing_2pt5 = outlet_staticpressure_blowing_2pt5(:,2);

array_size = size(distances_blowing_2pt5);
array_size = array_size(1);

integrand_array = zeros(array_size,1);

for n = 1:(array_size-1)

integrand_array(n) = 0.5*(distances_blowing_2pt5(n+1) - distances_blowing_2pt5(n))*(pressures_blowing_2pt5(n+1) + pressures_blowing_2pt5(n));

end

static_pressure_integrated_blowing_2pt5 = sum(integrand_array);


distances_blowing_3 = outlet_staticpressure_blowing_3(:,1);
pressures_blowing_3 = outlet_staticpressure_blowing_3(:,2);

array_size = size(distances_blowing_3);
array_size = array_size(1);

integrand_array = zeros(array_size,1);

for n = 1:(array_size-1)

integrand_array(n) = 0.5*(distances_blowing_3(n+1) - distances_blowing_3(n))*(pressures_blowing_3(n+1) + pressures_blowing_3(n));

end

static_pressure_integrated_blowing_3 = sum(integrand_array);


distances_blowing_3pt5 = outlet_staticpressure_blowing_3pt5(:,1);
pressures_blowing_3pt5 = outlet_staticpressure_blowing_3pt5(:,2);

array_size = size(distances_blowing_3pt5);
array_size = array_size(1);

integrand_array = zeros(array_size,1);

for n = 1:(array_size-1)

integrand_array(n) = 0.5*(distances_blowing_3pt5(n+1) - distances_blowing_3pt5(n))*(pressures_blowing_3pt5(n+1) + pressures_blowing_3pt5(n));

end

static_pressure_integrated_blowing_3pt5 = sum(integrand_array);


distances_blowing_4 = outlet_staticpressure_blowing_4(:,1);
pressures_blowing_4 = outlet_staticpressure_blowing_4(:,2);

array_size = size(distances_blowing_4);
array_size = array_size(1);

integrand_array = zeros(array_size,1);

for n = 1:(array_size-1)

integrand_array(n) = 0.5*(distances_blowing_4(n+1) - distances_blowing_4(n))*(pressures_blowing_4(n+1) + pressures_blowing_4(n));

end

static_pressure_integrated_blowing_4 = sum(integrand_array);



%%%%%%%%%%%LOSS VS TE BLOWING - COMBINING INTO ONE ARRAY%%%%%%%%%%%%%%%%%%%


total_pressures_integrated_blowing(1) = total_pressure_integrated_blowing_0;
total_pressures_integrated_blowing(2) = total_pressure_integrated_blowing_0pt5;
total_pressures_integrated_blowing(3) = total_pressure_integrated_blowing_1;
total_pressures_integrated_blowing(4) = total_pressure_integrated_blowing_1pt5;
total_pressures_integrated_blowing(5) = total_pressure_integrated_blowing_2;
total_pressures_integrated_blowing(6) = total_pressure_integrated_blowing_2pt5;
total_pressures_integrated_blowing(7) = total_pressure_integrated_blowing_3;
total_pressures_integrated_blowing(8) = total_pressure_integrated_blowing_3pt5;
total_pressures_integrated_blowing(9) = total_pressure_integrated_blowing_4;

static_pressures_integrated_blowing(1) = static_pressure_integrated_blowing_0;
static_pressures_integrated_blowing(2) = static_pressure_integrated_blowing_0pt5;
static_pressures_integrated_blowing(3) = static_pressure_integrated_blowing_1;
static_pressures_integrated_blowing(4) = static_pressure_integrated_blowing_1pt5;
static_pressures_integrated_blowing(5) = static_pressure_integrated_blowing_2;
static_pressures_integrated_blowing(6) = static_pressure_integrated_blowing_2pt5;
static_pressures_integrated_blowing(7) = static_pressure_integrated_blowing_3;
static_pressures_integrated_blowing(8) = static_pressure_integrated_blowing_3pt5;
static_pressures_integrated_blowing(9) = static_pressure_integrated_blowing_4;

loss_array_size = size(total_pressures_integrated_blowing);

loss_array_blowing = 0.01*100*(ones(loss_array_size) - (total_pressures_integrated_blowing)./(4325590));
%4325590*0.0588106 is the upstream TP times the inlet width

energy_loss_array_blowing = (0.01*100)*(ones(loss_array_size) - ((ones(loss_array_size) - ((static_pressures_integrated_blowing)./(total_pressures_integrated_blowing)).^0.285714)./(ones(loss_array_size) - ((static_pressures_integrated_blowing)./(4325590)).^0.285714)));
%needs a fudge factor of 0.01 for some reason







%%%%%%%%%%%%%%%%%%%%%%%%TURNING ANGLES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

turning_angles_00 = -360/(2*pi)*atan(yvel00(:,2)./xvel00(:,2));
turning_angles_01 = -360/(2*pi)*atan(yvel01(:,2)./xvel01(:,2));
turning_angles_02 = -360/(2*pi)*atan(yvel02(:,2)./xvel02(:,2));
turning_angles_03 = -360/(2*pi)*atan(yvel03(:,2)./xvel03(:,2));
turning_angles_04 = -360/(2*pi)*atan(yvel04(:,2)./xvel04(:,2));
turning_angles_05 = -360/(2*pi)*atan(yvel05(:,2)./xvel05(:,2));
turning_angles_06 = -360/(2*pi)*atan(yvel06(:,2)./xvel06(:,2));
turning_angles_07 = -360/(2*pi)*atan(yvel07(:,2)./xvel07(:,2));
turning_angles_08 = -360/(2*pi)*atan(yvel08(:,2)./xvel08(:,2));
turning_angles_09 = -360/(2*pi)*atan(yvel09(:,2)./xvel09(:,2));
turning_angles_10 = -360/(2*pi)*atan(yvel10(:,2)./xvel10(:,2));

turning_angle_means = zeros(11,1);
turning_angle_means(1) = mean(turning_angles_00);
turning_angle_means(2) = mean(turning_angles_01);
turning_angle_means(3) = mean(turning_angles_02);
turning_angle_means(4) = mean(turning_angles_03);
turning_angle_means(5) = mean(turning_angles_04);
turning_angle_means(6) = mean(turning_angles_05);
turning_angle_means(7) = mean(turning_angles_06);
turning_angle_means(8) = mean(turning_angles_07);
turning_angle_means(9) = mean(turning_angles_08);
turning_angle_means(10) = mean(turning_angles_09);
turning_angle_means(11) = mean(turning_angles_10);


turning_angles_ps_0 = -360/(2*pi)*atan(yvel00(:,2)./xvelps0(:,2));
turning_angles_ps_1 = -360/(2*pi)*atan(yvel01(:,2)./xvelps1(:,2));
turning_angles_ps_2 = -360/(2*pi)*atan(yvel02(:,2)./xvelps2(:,2));
turning_angles_ps_3 = -360/(2*pi)*atan(yvel03(:,2)./xvelps3(:,2));
turning_angles_ps_4 = -360/(2*pi)*atan(yvel04(:,2)./xvelps4(:,2));

turning_angle_ps_means = zeros(5,1);
turning_angle_ps_means(1) = mean(turning_angles_ps_0);
turning_angle_ps_means(2) = mean(turning_angles_ps_1);
turning_angle_ps_means(3) = mean(turning_angles_ps_2);
turning_angle_ps_means(4) = mean(turning_angles_ps_3);
turning_angle_ps_means(5) = mean(turning_angles_ps_4);





%%%%%%%%%%%%%%%%%%%%%%%%%LOSSES%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%4325590 is the absolute upstream total pressure
%2416530.7 is the absolute outlet static pressure at design
%4285900 is the operating pressure of the Fluent simulation, and is needed
%to offset the values of outlet static pressure reported by Fluent's gauge

losses_tp_00 = (4325590 - (tp00(:,2)+4285900) )./(4325590 - 2416530.7);
losses_tp_01 = (4325590 - (tp01(:,2)+4285900) )./(4325590 - 2416530.7);
losses_tp_02 = (4325590 - (tp02(:,2)+4285900) )./(4325590 - 2416530.7);
losses_tp_03 = (4325590 - (tp03(:,2)+4285900) )./(4325590 - 2416530.7);
losses_tp_04 = (4325590 - (tp04(:,2)+4285900) )./(4325590 - 2416530.7);
losses_tp_05 = (4325590 - (tp05(:,2)+4285900) )./(4325590 - 2416530.7);
losses_tp_06 = (4325590 - (tp06(:,2)+4285900) )./(4325590 - 2416530.7);
losses_tp_07 = (4325590 - (tp07(:,2)+4285900) )./(4325590 - 2416530.7);
losses_tp_08 = (4325590 - (tp08(:,2)+4285900) )./(4325590 - 2416530.7);
losses_tp_09 = (4325590 - (tp09(:,2)+4285900) )./(4325590 - 2416530.7);
losses_tp_10 = (4325590 - (tp10(:,2)+4285900) )./(4325590 - 2416530.7);

losses_tp_means = zeros(11,1);
losses_tp_means(1) = mean(losses_tp_00);
losses_tp_means(2) = mean(losses_tp_01);
losses_tp_means(3) = mean(losses_tp_02);
losses_tp_means(4) = mean(losses_tp_03);
losses_tp_means(5) = mean(losses_tp_04);
losses_tp_means(6) = mean(losses_tp_05);
losses_tp_means(7) = mean(losses_tp_06);
losses_tp_means(8) = mean(losses_tp_07);
losses_tp_means(9) = mean(losses_tp_08);
losses_tp_means(10) = mean(losses_tp_09);
losses_tp_means(11) = mean(losses_tp_10);

losses_ke_00 = ( (4325590./(tp00(:,2)+4285900)).^0.2857 - 1 )/( (4325590/2416530.7)^0.2857 - 1 );
losses_ke_01 = ( (4325590./(tp01(:,2)+4285900)).^0.2857 - 1 )/( (4325590/2416530.7)^0.2857 - 1 );
losses_ke_02 = ( (4325590./(tp02(:,2)+4285900)).^0.2857 - 1 )/( (4325590/2416530.7)^0.2857 - 1 );
losses_ke_03 = ( (4325590./(tp03(:,2)+4285900)).^0.2857 - 1 )/( (4325590/2416530.7)^0.2857 - 1 );
losses_ke_04 = ( (4325590./(tp04(:,2)+4285900)).^0.2857 - 1 )/( (4325590/2416530.7)^0.2857 - 1 );
losses_ke_05 = ( (4325590./(tp05(:,2)+4285900)).^0.2857 - 1 )/( (4325590/2416530.7)^0.2857 - 1 );
losses_ke_06 = ( (4325590./(tp06(:,2)+4285900)).^0.2857 - 1 )/( (4325590/2416530.7)^0.2857 - 1 );
losses_ke_07 = ( (4325590./(tp07(:,2)+4285900)).^0.2857 - 1 )/( (4325590/2416530.7)^0.2857 - 1 );
losses_ke_08 = ( (4325590./(tp08(:,2)+4285900)).^0.2857 - 1 )/( (4325590/2416530.7)^0.2857 - 1 );
losses_ke_09 = ( (4325590./(tp09(:,2)+4285900)).^0.2857 - 1 )/( (4325590/2416530.7)^0.2857 - 1 );
losses_ke_10 = ( (4325590./(tp10(:,2)+4285900)).^0.2857 - 1 )/( (4325590/2416530.7)^0.2857 - 1 );

losses_ke_means = zeros(11,1);
losses_ke_means(1) = mean(losses_ke_00);
losses_ke_means(2) = mean(losses_ke_01);
losses_ke_means(3) = mean(losses_ke_02);
losses_ke_means(4) = mean(losses_ke_03);
losses_ke_means(5) = mean(losses_ke_04);
losses_ke_means(6) = mean(losses_ke_05);
losses_ke_means(7) = mean(losses_ke_06);
losses_ke_means(8) = mean(losses_ke_07);
losses_ke_means(9) = mean(losses_ke_08);
losses_ke_means(10) = mean(losses_ke_09);
losses_ke_means(11) = mean(losses_ke_10);


losses_tp_ps_0 = (4325590 - (tpps0(:,2)+4285900) )./(4325590 - 2416530.7);
losses_tp_ps_1 = (4325590 - (tpps1(:,2)+4285900) )./(4325590 - 2416530.7);
losses_tp_ps_2 = (4325590 - (tpps2(:,2)+4285900) )./(4325590 - 2416530.7);
losses_tp_ps_3 = (4325590 - (tpps3(:,2)+4285900) )./(4325590 - 2416530.7);
losses_tp_ps_4 = (4325590 - (tpps4(:,2)+4285900) )./(4325590 - 2416530.7);

losses_tp_ps_means = zeros(5,1);
losses_tp_ps_means(1) = mean(losses_tp_ps_0);
losses_tp_ps_means(2) = mean(losses_tp_ps_1);
losses_tp_ps_means(3) = mean(losses_tp_ps_2);
losses_tp_ps_means(4) = mean(losses_tp_ps_3);
losses_tp_ps_means(5) = mean(losses_tp_ps_4);

losses_ke_ps_0 = ( (4325590./(tpps0(:,2)+4285900)).^0.2857 - 1 )/( (4325590/2416530.7)^0.2857 - 1 );
losses_ke_ps_1 = ( (4325590./(tpps1(:,2)+4285900)).^0.2857 - 1 )/( (4325590/2416530.7)^0.2857 - 1 );
losses_ke_ps_2 = ( (4325590./(tpps2(:,2)+4285900)).^0.2857 - 1 )/( (4325590/2416530.7)^0.2857 - 1 );
losses_ke_ps_3 = ( (4325590./(tpps3(:,2)+4285900)).^0.2857 - 1 )/( (4325590/2416530.7)^0.2857 - 1 );
losses_ke_ps_4 = ( (4325590./(tpps4(:,2)+4285900)).^0.2857 - 1 )/( (4325590/2416530.7)^0.2857 - 1 );

losses_ke_ps_means = zeros(5,1);
losses_ke_ps_means(1) = mean(losses_ke_ps_0);
losses_ke_ps_means(2) = mean(losses_ke_ps_1);
losses_ke_ps_means(3) = mean(losses_ke_ps_2);
losses_ke_ps_means(4) = mean(losses_ke_ps_3);
losses_ke_ps_means(5) = mean(losses_ke_ps_4);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%PLOTS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Default figure settings
figure_width = 4;     % Width in inches
figure_height = 4;    % Height in inches
axes_line_width = 0.75;    % AxesLineWidth
font_size = 14;      % Fontsize
line_width = 1.5;      % LineWidth
marker_size = 8;       % MarkerSize


%set(0, 'defaultfigureunits', 'centimeters')
%set(0, 'defaultfigureposition', [35 15 9 9])
%set(0, 'defaultaxesposition', [0.19 0.14 0.76 0.75])


figure(2)
%plot of normalised NGV capacity against NGV pressure ratio
plot(T900_capacities_2D(18:28,1), T900_capacities_2D(18:28,7)./0.000004781, 'k.')
hold on
plot(T900_capacities_2D(18:28,1), T900_capacities_2D(18:28,5)./0.000004781, 'r.')
plot(T900_capacities_2D(18:28,1), T900_capacities_2D(18:28,6)./0.000004781, 'm.')
plot(T900_capacities_2D(18:28,1), T900_capacities_2D(18:28,3)./0.000004781, 'g.')
plot(T900_capacities_2D(18:28,1), T900_capacities_2D(18:28,2)./0.000004781, 'b.')
plot(T900_capacities_2D(18:28,1), T900_capacities_2D(18:28,4)./0.000004781, 'c.')
plot(T900_capacities_2D_oversampled(:,1), T900_capacities_2D_oversampled(:,7)./0.000004781, 'k--')
plot(T900_capacities_2D_oversampled(:,1), T900_capacities_2D_oversampled(:,5)./0.000004781, 'r--')
plot(T900_capacities_2D_oversampled(:,1), T900_capacities_2D_oversampled(:,6)./0.000004781, 'm--')
plot(T900_capacities_2D_oversampled(:,1), T900_capacities_2D_oversampled(:,3)./0.000004781, 'g-')
plot(T900_capacities_2D_oversampled(:,1), T900_capacities_2D_oversampled(:,2)./0.000004781, 'b-')
plot(T900_capacities_2D_oversampled(:,1), T900_capacities_2D_oversampled(:,4)./0.000004781, 'c-')
xlim([1.5 2.5])
xlabel('NGV pressure ratio, $\frac{p_{01}}{p_2}$','Interpreter','latex')
ylabel('\Delta capacity, %')
design_line = xline(1.79,'-k',{'Design'},'FontName','Charter','FontSize',font_size);
design_line.LabelVerticalAlignment = 'bottom';
%legend('PNS04','PNN06', 'PNS03', 'KTA01', 'KSZ03', 'KVD04', 'Location', 'Best', 'FontName','Charter','FontSize',font_size)

%Here we set up the axes
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) figure_width*100, figure_height*100]); %<- Set size
set(gca, 'FontSize', font_size, 'LineWidth', axes_line_width); %<- Set properties
set(gca,'FontName','Charter','FontSize',font_size)

% Here we preserve the size of the image when we save it.
set(gcf,'InvertHardcopy','on');
set(gcf,'PaperUnits', 'inches');
papersize = get(gcf, 'PaperSize');
left = (papersize(1)- figure_width)/2;
bottom = (papersize(2)- figure_height)/2;
myfiguresize = [left, bottom, figure_width, figure_height];
set(gcf,'PaperPosition', myfiguresize);
print('../../figs/t900_2d_capacity_trends','-dpng','-r300');


figure(3)
%plot of Rolls-Royce 3D CFD data
plot(T900_capacities_3D(:,11), T900_capacities_3D(:,12)/15200.3/0.00001157, 'k.')
hold on
plot(T900_capacities_3D(:,7), T900_capacities_3D(:,8)/15200.3/0.00001157, 'r.')
plot(T900_capacities_3D(:,9), T900_capacities_3D(:,10)/15200.3/0.00001157, 'm.')
plot(T900_capacities_3D(:,3), T900_capacities_3D(:,4)/15200.3/0.00001157, 'g.')
plot(T900_capacities_3D(:,1), T900_capacities_3D(:,2)/15200.3/0.00001157, 'b.')
plot(T900_capacities_3D(:,5), T900_capacities_3D(:,6)/15200.3/0.00001157, 'c.')
plot(T900_capacity_3D_6_oversampled(:,1), T900_capacity_3D_6_oversampled(:,2)/15200.3/0.00001157, 'k--')
hold on
plot(T900_capacity_3D_4_oversampled(:,1), T900_capacity_3D_4_oversampled(:,2)/15200.3/0.00001157, 'r--')
plot(T900_capacity_3D_5_oversampled(:,1), T900_capacity_3D_5_oversampled(:,2)/15200.3/0.00001157, 'm--')
plot(T900_capacity_3D_2_oversampled(:,1), T900_capacity_3D_2_oversampled(:,2)/15200.3/0.00001157, 'g-')
plot(T900_capacity_3D_1_oversampled(:,1), T900_capacity_3D_1_oversampled(:,2)/15200.3/0.00001157, 'b-')
plot(T900_capacity_3D_3_oversampled(:,1), T900_capacity_3D_3_oversampled(:,2)/15200.3/0.00001157, 'c-')
xlim([1.5 2.5]) 
xlabel('NGV pressure ratio, $\frac{p_{01}}{p_2}$','Interpreter','latex')
ylabel('\Delta capacity, %')
design_line = xline(1.79,'-k',{'Design'},'FontName','Charter','FontSize',font_size);
design_line.LabelVerticalAlignment = 'bottom';
%legend('PNS04','PNN06', 'PNS03', 'KTA01', 'KSZ03', 'KVD04', 'Location', 'Best')

%Here we set up the axes
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) figure_width*100, figure_height*100]); %<- Set size
set(gca, 'FontSize', font_size, 'LineWidth', axes_line_width); %<- Set properties
set(gca,'FontName','Charter','FontSize',font_size)

% Here we preserve the size of the image when we save it.
set(gcf,'InvertHardcopy','on');
set(gcf,'PaperUnits', 'inches');
papersize = get(gcf, 'PaperSize');
left = (papersize(1)- figure_width)/2;
bottom = (papersize(2)- figure_height)/2;
myfiguresize = [left, bottom, figure_width, figure_height];
set(gcf,'PaperPosition', myfiguresize);
print('../../figs/t900_3d_capacity_trends','-dpng','-r300');


figure(4)
%2D Capacity versus throat width
%plot(100*throat_widths(1:3)./max(throat_widths) , 100*T900_capacities_2D_design(2:4)./max(T900_capacities_2D_design(2:7)), 'ro')
plot(100*throat_widths(1)./max(throat_widths) , 100*T900_capacities_2D_design(2)./max(T900_capacities_2D_design(2:7)), 'bo','MarkerFaceColor','blue','MarkerSize',10)
hold on
%plot(100*throat_widths(4:6)./max(throat_widths) , 100*T900_capacities_2D_design(5:7)./max(T900_capacities_2D_design(2:7)), 'bo')
plot(100*throat_widths(2)./max(throat_widths) , 100*T900_capacities_2D_design(3)./max(T900_capacities_2D_design(2:7)), 'go','MarkerFaceColor','green','MarkerSize',10)
plot(100*throat_widths(3)./max(throat_widths) , 100*T900_capacities_2D_design(4)./max(T900_capacities_2D_design(2:7)), 'co','MarkerFaceColor','cyan','MarkerSize',10)
plot(100*throat_widths(4)./max(throat_widths) , 100*T900_capacities_2D_design(5)./max(T900_capacities_2D_design(2:7)), 'ro','MarkerFaceColor','red','MarkerSize',10)
plot(100*throat_widths(5)./max(throat_widths) , 100*T900_capacities_2D_design(6)./max(T900_capacities_2D_design(2:7)), 'mo','MarkerFaceColor','magenta','MarkerSize',10)
plot(100*throat_widths(6)./max(throat_widths) , 100*T900_capacities_2D_design(7)./max(T900_capacities_2D_design(2:7)), 'ko','MarkerFaceColor','black','MarkerSize',10)
plot(100*throat_width_Mskew_mean/max(throat_widths), 100*T900_capacity_Mskew_mean/max(T900_capacities_2D_design(2:7)), 'k*')
plot(100*throat_width_EP1_mean/max(throat_widths), 100*T900_capacity_EP1_mean/max(T900_capacities_2D_design(2:7)), 'k*')
plot(1.4*gradient_1(:,1) + 100*throat_width_Mskew_mean/max(throat_widths), 1.4*gradient_1(:,2) + 100*T900_capacity_Mskew_mean/max(T900_capacities_2D_design(2:7)), 'k-')
plot(1000*gradient_1(:,1) + 100*throat_width_EP1_mean/max(throat_widths), 1000*gradient_1(:,2) + 100*T900_capacity_EP1_mean/max(T900_capacities_2D_design(2:7)), 'k')
plot(-1000*gradient_1(:,1) + 100*throat_width_Mskew_mean/max(throat_widths), -1000*gradient_1(:,2) + 100*T900_capacity_Mskew_mean/max(T900_capacities_2D_design(2:7)), 'k-')
plot(-0.8*gradient_1(:,1) + 100*throat_width_EP1_mean/max(throat_widths), -0.8*gradient_1(:,2) + 100*T900_capacity_EP1_mean/max(T900_capacities_2D_design(2:7)), 'k')
axis([100*0.99*min(throat_widths)/max(throat_widths) 100*1.01 100*0.99*min(T900_capacities_2D_design(2:7))/max(T900_capacities_2D_design(2:7)) 100*1.01])
xlabel('\Delta throat width, %')
ylabel('\Delta capacity, %')
%legend('M-skewed', 'EP1', 'Avg. M-skewed vane', 'Avg. EP1 vane', 'Line of gradient 1', 'Line of gradient 1', 'Location', 'NorthWest')
%text(100*throat_widths(:)./max(throat_widths) + 100*0.002, 100*T900_capacities_2D_design(2:7)./max(T900_capacities_2D_design(2:7)), {'KSZ03', 'KTA01', 'KVD04', 'PNN06', 'PNS03', 'PNS04'},'FontName','Charter','FontSize',font_size)
%text(100*throat_widths(1)./max(throat_widths) + 100*0.002, 100*T900_capacities_2D_design(2)./max(T900_capacities_2D_design(2:7)),'KSZ03','FontName','Charter','FontSize',font_size,'VerticalAlignment','top','HorizontalAlignment','left')
%text(100*throat_widths(2)./max(throat_widths) + 100*0.002, 100*T900_capacities_2D_design(3)./max(T900_capacities_2D_design(2:7)),'KTA01','FontName','Charter','FontSize',font_size,'VerticalAlignment','bottom','HorizontalAlignment','right')
%text(100*throat_widths(3)./max(throat_widths) + 100*0.002, 100*T900_capacities_2D_design(4)./max(T900_capacities_2D_design(2:7)),'KVD04','FontName','Charter','FontSize',font_size,'VerticalAlignment','top','HorizontalAlignment','left')
%text(100*throat_widths(4)./max(throat_widths) + 100*0.002, 100*T900_capacities_2D_design(5)./max(T900_capacities_2D_design(2:7)),'PNN06','FontName','Charter','FontSize',font_size,'VerticalAlignment','top','HorizontalAlignment','left')
%text(100*throat_widths(5)./max(throat_widths) + 100*0.002, 100*T900_capacities_2D_design(6)./max(T900_capacities_2D_design(2:7)),'PNS03','FontName','Charter','FontSize',font_size,'VerticalAlignment','bottom','HorizontalAlignment','right')
%text(100*throat_widths(6)./max(throat_widths) + 100*0.002, 100*T900_capacities_2D_design(7)./max(T900_capacities_2D_design(2:7)),'PNS04','FontName','Charter','FontSize',font_size,'VerticalAlignment','top','HorizontalAlignment','left')
%Here we set up the axes
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) figure_width*100, figure_height*100]); %<- Set size
set(gca, 'FontSize', font_size, 'LineWidth', axes_line_width); %<- Set properties
set(gca,'FontName','Charter','FontSize',font_size)
% Here we preserve the size of the image when we save it.
set(gcf,'InvertHardcopy','on');
set(gcf,'PaperUnits', 'inches');
papersize = get(gcf, 'PaperSize');
left = (papersize(1)- figure_width)/2;
bottom = (papersize(2)- figure_height)/2;
myfiguresize = [left, bottom, figure_width, figure_height];
set(gcf,'PaperPosition', myfiguresize);
print('../../figs/t900_2d_capacities_vs_throat_widths','-dpng','-r300');


figure(5)
%3D Capacity versus throat area
%plot(100*throat_areas_digitised(1:3)./max(throat_areas_digitised) , 100*T900_capacities_3D_design(2:4)./max(T900_capacities_3D_design(2:7)), 'ro')
plot(100*throat_areas_digitised(1)./max(throat_areas_digitised) , 100*T900_capacities_3D_design(2)./max(T900_capacities_3D_design(2:7)), 'bo','MarkerFaceColor','blue','MarkerSize',10)
hold on
plot(100*throat_areas_digitised(4:6)./max(throat_areas_digitised) , 100*T900_capacities_3D_design(5:7)./max(T900_capacities_3D_design(2:7)), 'bo')
plot(100*throat_areas_digitised(2)./max(throat_areas_digitised) , 100*T900_capacities_3D_design(3)./max(T900_capacities_3D_design(2:7)), 'go','MarkerFaceColor','green','MarkerSize',10)
plot(100*throat_areas_digitised(3)./max(throat_areas_digitised) , 100*T900_capacities_3D_design(4)./max(T900_capacities_3D_design(2:7)), 'co','MarkerFaceColor','cyan','MarkerSize',10)
plot(100*throat_areas_digitised(4)./max(throat_areas_digitised) , 100*T900_capacities_3D_design(5)./max(T900_capacities_3D_design(2:7)), 'ro','MarkerFaceColor','red','MarkerSize',10)
plot(100*throat_areas_digitised(5)./max(throat_areas_digitised) , 100*T900_capacities_3D_design(6)./max(T900_capacities_3D_design(2:7)), 'mo','MarkerFaceColor','magenta','MarkerSize',10)
plot(100*throat_areas_digitised(6)./max(throat_areas_digitised) , 100*T900_capacities_3D_design(7)./max(T900_capacities_3D_design(2:7)), 'ko','MarkerFaceColor','black','MarkerSize',10)
plot(100*throat_width_Mskew_mean_3D/max(throat_areas_digitised), 100*T900_capacity_Mskew_mean_3D/max(T900_capacities_3D_design(2:7)), 'k*')
plot(100*throat_width_EP1_mean_3D/max(throat_areas_digitised), 100*T900_capacity_EP1_mean_3D/max(T900_capacities_3D_design(2:7)), 'k*')
plot(1.4*gradient_1(:,1) + 100*throat_width_Mskew_mean_3D/max(throat_areas_digitised), 1.4*gradient_1(:,2) + 100*T900_capacity_Mskew_mean_3D/max(T900_capacities_3D_design(2:7)), 'k-')
plot(1000*gradient_1(:,1) + 100*throat_width_EP1_mean_3D/max(throat_areas_digitised), 1000*gradient_1(:,2) + 100*T900_capacity_EP1_mean_3D/max(T900_capacities_3D_design(2:7)), 'k')
plot(-1000*gradient_1(:,1) + 100*throat_width_Mskew_mean_3D/max(throat_areas_digitised), -1000*gradient_1(:,2) + 100*T900_capacity_Mskew_mean_3D/max(T900_capacities_3D_design(2:7)), 'k-')
plot(-0.8*gradient_1(:,1) + 100*throat_width_EP1_mean_3D/max(throat_areas_digitised), -0.8*gradient_1(:,2) + 100*T900_capacity_EP1_mean_3D/max(T900_capacities_3D_design(2:7)), 'k')
axis([100*0.99*min(throat_areas_digitised)/max(throat_areas_digitised) 100*1.01 100*0.99*min(T900_capacities_3D_design(2:7))/max(T900_capacities_3D_design(2:7)) 100*1.01])
xlabel('\Delta throat area, %')
ylabel('\Delta capacity, %')
%legend('M-skewed', 'EP1', 'Avg. M-skewed vane', 'Avg. EP1 vane', 'Line of gradient 1', 'Line of gradient 1', 'Location', 'NorthWest')
%text(100*throat_areas_digitised(:)./max(throat_areas_digitised) + 100*0.002, 100*T900_capacities_3D_design(2:7)./max(T900_capacities_3D_design(2:7)), {'KSZ03', 'KTA01', 'KVD04', 'PNN06', 'PNS03', 'PNS04'},'FontName','Charter','FontSize',font_size)
%text(100*throat_areas_digitised(1)./max(throat_areas_digitised) + 100*0.002, 100*T900_capacities_3D_design(2)./max(T900_capacities_3D_design(2:7)),'KSZ03','FontName','Charter','FontSize',font_size,'VerticalAlignment','top','HorizontalAlignment','left')
%text(100*throat_areas_digitised(2)./max(throat_areas_digitised) + 100*0.002, 100*T900_capacities_3D_design(3)./max(T900_capacities_3D_design(2:7)),'KTA01','FontName','Charter','FontSize',font_size,'VerticalAlignment','bottom','HorizontalAlignment','right')
%text(100*throat_areas_digitised(3)./max(throat_areas_digitised) + 100*0.002, 100*T900_capacities_3D_design(4)./max(T900_capacities_3D_design(2:7)),'KVD04','FontName','Charter','FontSize',font_size,'VerticalAlignment','top','HorizontalAlignment','left')
%text(100*throat_areas_digitised(4)./max(throat_areas_digitised) + 100*0.002, 100*T900_capacities_3D_design(5)./max(T900_capacities_3D_design(2:7)),'PNN06','FontName','Charter','FontSize',font_size,'VerticalAlignment','top','HorizontalAlignment','left')
%text(100*throat_areas_digitised(5)./max(throat_areas_digitised) + 100*0.002, 100*T900_capacities_3D_design(6)./max(T900_capacities_3D_design(2:7)),'PNS03','FontName','Charter','FontSize',font_size,'VerticalAlignment','bottom','HorizontalAlignment','right')
%text(100*throat_areas_digitised(6)./max(throat_areas_digitised) + 100*0.002, 100*T900_capacities_3D_design(7)./max(T900_capacities_3D_design(2:7)),'PNS04','FontName','Charter','FontSize',font_size,'VerticalAlignment','top','HorizontalAlignment','left')
%Here we set up the axes
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) figure_width*100, figure_height*100]); %<- Set size
set(gca, 'FontSize', font_size, 'LineWidth', axes_line_width); %<- Set properties
set(gca,'FontName','Charter','FontSize',font_size)
% Here we preserve the size of the image when we save it.
set(gcf,'InvertHardcopy','on');
set(gcf,'PaperUnits', 'inches');
papersize = get(gcf, 'PaperSize');
left = (papersize(1)- figure_width)/2;
bottom = (papersize(2)- figure_height)/2;
myfiguresize = [left, bottom, figure_width, figure_height];
set(gcf,'PaperPosition', myfiguresize);
print('../../figs/t900_3D_capacities_vs_throat_areas','-dpng','-r300');

throat_widths_effective = [0.0128,0.013,0.0127,0.0133,0.0133,0.0132]';
throat_widths_effective_mean = mean(throat_widths_effective);
T900_capacities_2D_design_mean = mean(T900_capacities_2D_design);

gradient = [0 0 ; 110 110];


figure(6)
%2D Capacity versus effective throat width

plot(100*throat_widths_effective(1)./max(throat_widths_effective) , 100*T900_capacities_2D_design(2)./max(T900_capacities_2D_design(2:7)), 'bo','MarkerFaceColor','blue','MarkerSize',10)
hold on
plot(100*throat_widths_effective(2)./max(throat_widths_effective) , 100*T900_capacities_2D_design(3)./max(T900_capacities_2D_design(2:7)), 'go','MarkerFaceColor','green','MarkerSize',10)
plot(100*throat_widths_effective(3)./max(throat_widths_effective) , 100*T900_capacities_2D_design(4)./max(T900_capacities_2D_design(2:7)), 'co','MarkerFaceColor','cyan','MarkerSize',10)
plot(100*throat_widths_effective(4)./max(throat_widths_effective) , 100*T900_capacities_2D_design(5)./max(T900_capacities_2D_design(2:7)), 'ro','MarkerFaceColor','red','MarkerSize',10)
plot(100*throat_widths_effective(5)./max(throat_widths_effective) , 100*T900_capacities_2D_design(6)./max(T900_capacities_2D_design(2:7)), 'mo','MarkerFaceColor','magenta','MarkerSize',10)
plot(100*throat_widths_effective(6)./max(throat_widths_effective) , 100*T900_capacities_2D_design(7)./max(T900_capacities_2D_design(2:7)), 'ko','MarkerFaceColor','black','MarkerSize',10)
plot(100*throat_widths_effective_mean/max(throat_widths_effective), 100*T900_capacities_2D_design_mean/max(T900_capacities_2D_design(2:7)), 'k*')
%plot(1.4*gradient_1(:,1) + 100*throat_widths_effective_mean/max(throat_widths_effective), 1.4*gradient_1(:,2) + 100*T900_capacities_2D_design_mean/max(T900_capacities_2D_design(2:7)), 'r-')
%plot(1000*gradient_1(:,1) + 100*throat_widths_effective_mean/max(throat_widths_effective), 1000*gradient_1(:,2) + 100*T900_capacities_2D_design_mean/max(T900_capacities_2D_design(2:7)), 'g')
%plot(-1000*gradient_1(:,1) + 100*throat_widths_effective_mean/max(throat_widths_effective), -1000*gradient_1(:,2) + 100*T900_capacities_2D_design_mean/max(T900_capacities_2D_design(2:7)), 'b-')
%plot(-0.8*gradient_1(:,1) + 100*throat_widths_effective_mean/max(throat_widths_effective), -0.8*gradient_1(:,2) + 100*T900_capacities_2D_design_mean/max(T900_capacities_2D_design(2:7)), 'c')
%plot(gradient(:,1) + 100*throat_widths_effective_mean/max(throat_widths_effective), -gradient(:,2) + 100*T900_capacities_2D_design_mean/max(T900_capacities_2D_design(2:7)), 'c')
plot(gradient(:,1),gradient(:,2),'k-')

axis([100*0.99*min(throat_widths_effective)/max(throat_widths_effective) 100*1.01 100*0.99*min(T900_capacities_2D_design(2:7))/max(T900_capacities_2D_design(2:7)) 100*1.01])
xlabel('\Delta effective throat width, %')
ylabel('\Delta capacity, %')
%legend('M-skewed', 'EP1', 'Avg. M-skewed vane', 'Avg. EP1 vane', 'Line of gradient 1', 'Line of gradient 1', 'Location', 'NorthWest')
%text(100*throat_widths(:)./max(throat_widths) + 100*0.002, 100*T900_capacities_2D_design(2:7)./max(T900_capacities_2D_design(2:7)), {'KSZ03', 'KTA01', 'KVD04', 'PNN06', 'PNS03', 'PNS04'},'FontName','Charter','FontSize',font_size)
%text(100*throat_widths(1)./max(throat_widths) + 100*0.002, 100*T900_capacities_2D_design(2)./max(T900_capacities_2D_design(2:7)),'KSZ03','FontName','Charter','FontSize',font_size,'VerticalAlignment','top','HorizontalAlignment','left')
%text(100*throat_widths(2)./max(throat_widths) + 100*0.002, 100*T900_capacities_2D_design(3)./max(T900_capacities_2D_design(2:7)),'KTA01','FontName','Charter','FontSize',font_size,'VerticalAlignment','bottom','HorizontalAlignment','right')
%text(100*throat_widths(3)./max(throat_widths) + 100*0.002, 100*T900_capacities_2D_design(4)./max(T900_capacities_2D_design(2:7)),'KVD04','FontName','Charter','FontSize',font_size,'VerticalAlignment','top','HorizontalAlignment','left')
%text(100*throat_widths(4)./max(throat_widths) + 100*0.002, 100*T900_capacities_2D_design(5)./max(T900_capacities_2D_design(2:7)),'PNN06','FontName','Charter','FontSize',font_size,'VerticalAlignment','top','HorizontalAlignment','left')
%text(100*throat_widths(5)./max(throat_widths) + 100*0.002, 100*T900_capacities_2D_design(6)./max(T900_capacities_2D_design(2:7)),'PNS03','FontName','Charter','FontSize',font_size,'VerticalAlignment','bottom','HorizontalAlignment','right')
%text(100*throat_widths(6)./max(throat_widths) + 100*0.002, 100*T900_capacities_2D_design(7)./max(T900_capacities_2D_design(2:7)),'PNS04','FontName','Charter','FontSize',font_size,'VerticalAlignment','top','HorizontalAlignment','left')

%Here we set up the axes
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) figure_width*100, figure_height*100]); %<- Set size
set(gca, 'FontSize', font_size, 'LineWidth', axes_line_width); %<- Set properties
set(gca,'FontName','Charter','FontSize',font_size)
% Here we preserve the size of the image when we save it.
set(gcf,'InvertHardcopy','on');
set(gcf,'PaperUnits', 'inches');
papersize = get(gcf, 'PaperSize');
left = (papersize(1)- figure_width)/2;
bottom = (papersize(2)- figure_height)/2;
myfiguresize = [left, bottom, figure_width, figure_height];
set(gcf,'PaperPosition', myfiguresize);
print('../../figs/t900_2d_capacities_vs_effective_throat_widths','-dpng','-r300');


figure(7)
%realistic TEs with SS cutbacks

plot(T900_TE_capacities_2D(:,1), 100*T900_TE_capacities_2D(:,2)./T900_TE_capacities_2D(18,2), 'k.')
hold on
%plot(T900_TE_capacities_2D(:,1), 100*T900_TE_capacities_2D(:,3)./T900_TE_capacities_2D(18,2), 'kx')
plot(T900_TE_capacities_2D(:,1), 100*T900_TE_capacities_2D(:,4)./T900_TE_capacities_2D(18,2), 'r.')
%plot(T900_TE_capacities_2D(:,1), 100*T900_TE_capacities_2D(:,5)./T900_TE_capacities_2D(18,2), 'kh')
plot(T900_TE_capacities_2D(:,1), 100*T900_TE_capacities_2D(:,6)./T900_TE_capacities_2D(18,2), 'g.')
%plot(T900_TE_capacities_2D(:,1), 100*T900_TE_capacities_2D(:,7)./T900_TE_capacities_2D(18,2), 'ko')
plot(T900_TE_capacities_2D(:,1), 100*T900_TE_capacities_2D(:,8)./T900_TE_capacities_2D(18,2), 'b.')
%plot(T900_TE_capacities_2D(:,1), 100*T900_TE_capacities_2D(:,9)./T900_TE_capacities_2D(18,2), 'kp')
plot(T900_TE_capacities_2D(:,1), 100*T900_TE_capacities_2D(:,10)./T900_TE_capacities_2D(18,2), 'c.')
%plot(T900_TE_capacities_2D(:,1), 100*T900_TE_capacities_2D(:,11)./T900_TE_capacities_2D(18,2), 'kd')
plot(T900_TE_capacities_2D(:,1), 100*T900_TE_capacities_2D(:,12)./T900_TE_capacities_2D(18,2), 'm.')

plot(T900_TE_capacities_2D_oversampled(:,1), 100*T900_TE_capacities_2D_oversampled(:,2)./T900_TE_capacities_2D(18,2), 'k-')
%plot(T900_TE_capacities_2D_oversampled(:,1), 100*T900_TE_capacities_2D_oversampled(:,3)./T900_TE_capacities_2D(18,2), 'k-.')
plot(T900_TE_capacities_2D_oversampled(:,1), 100*T900_TE_capacities_2D_oversampled(:,4)./T900_TE_capacities_2D(18,2), 'r-')
%plot(T900_TE_capacities_2D_oversampled(:,1), 100*T900_TE_capacities_2D_oversampled(:,5)./T900_TE_capacities_2D(18,2), 'k-.')
plot(T900_TE_capacities_2D_oversampled(:,1), 100*T900_TE_capacities_2D_oversampled(:,6)./T900_TE_capacities_2D(18,2), 'g-')
%plot(T900_TE_capacities_2D_oversampled(:,1), 100*T900_TE_capacities_2D_oversampled(:,7)./T900_TE_capacities_2D(18,2), 'k-.')
plot(T900_TE_capacities_2D_oversampled(:,1), 100*T900_TE_capacities_2D_oversampled(:,8)./T900_TE_capacities_2D(18,2), 'b-')
%plot(T900_TE_capacities_2D_oversampled(:,1), 100*T900_TE_capacities_2D_oversampled(:,9)./T900_TE_capacities_2D(18,2), 'k-.')
plot(T900_TE_capacities_2D_oversampled(:,1), 100*T900_TE_capacities_2D_oversampled(:,10)./T900_TE_capacities_2D(18,2), 'c-')
%plot(T900_TE_capacities_2D_oversampled(:,1), 100*T900_TE_capacities_2D_oversampled(:,11)./T900_TE_capacities_2D(18,2), 'k-.')
plot(T900_TE_capacities_2D_oversampled(:,1), 100*T900_TE_capacities_2D_oversampled(:,12)./T900_TE_capacities_2D(18,2), 'm-')

xlim([1.5 2.5])
xlabel('NGV pressure ratio, $\frac{p_{01}}{p_2}$','Interpreter','latex')
ylabel('\Delta capacity, %')
design_line = xline(1.79,'-k',{'Design'},'FontName','Charter','FontSize',font_size);
design_line.LabelVerticalAlignment = 'bottom';
legend('0%', '20%', '40%', '60%', '80%', '100%', 'Location', 'southeast')

%Here we set up the axes
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) figure_width*100, figure_height*100]); %<- Set size
set(gca, 'FontSize', font_size, 'LineWidth', axes_line_width); %<- Set properties
set(gca,'FontName','Charter','FontSize',font_size)
% Here we preserve the size of the image when we save it.
set(gcf,'InvertHardcopy','on');
set(gcf,'PaperUnits', 'inches');
papersize = get(gcf, 'PaperSize');
left = (papersize(1)- figure_width)/2;
bottom = (papersize(2)- figure_height)/2;
myfiguresize = [left, bottom, figure_width, figure_height];
set(gcf,'PaperPosition', myfiguresize);
print('../../figs/ss_cutbacks_vs_capacities_trends','-dpng','-r300');


figure(8)
%capacity vs. cutback at various PRs

plot(0:10:100, 100*T900_TE_capacities_2D_oversampled(2105,2:12)./T900_TE_capacities_2D_oversampled(2105,2)-100, 'ro-')
hold on
%plot(0:10:100, 100*T900_TE_capacities_2D_oversampled(6670,2:12)./T900_TE_capacities_2D_oversampled(3354,2), 'kx-')
plot(0:10:100, 100*T900_TE_capacities_2D_oversampled(3354,2:12)./T900_TE_capacities_2D_oversampled(3354,2)-100, 'bo-')
plot(0:10:100, 100*T900_TE_capacities_2D_oversampled(10001,2:12)./T900_TE_capacities_2D_oversampled(10001,2)-100, 'mo-')
%plot(0:10:100, 100*T900_TE_capacities_2D_oversampled(37,2:12)./T900_TE_capacities_2D_oversampled(3354,2), 'k^-')

xlabel('Cutback amount, %')
ylabel('\Delta capacity, %')
legend('PR = 1.50', 'PR = 1.79 (design)', 'PR = 3.33', 'Position',[0.64 0.28 0.05 0.14])

%Here we set up the axes
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) figure_width*100, figure_height*100]); %<- Set size
set(gca, 'FontSize', font_size, 'LineWidth', axes_line_width); %<- Set properties
set(gca,'FontName','Charter','FontSize',font_size)
% Here we preserve the size of the image when we save it.
set(gcf,'InvertHardcopy','on');
set(gcf,'PaperUnits', 'inches');
papersize = get(gcf, 'PaperSize');
left = (papersize(1)- figure_width)/2;
bottom = (papersize(2)- figure_height)/2;
myfiguresize = [left, bottom, figure_width, figure_height];
set(gcf,'PaperPosition', myfiguresize);
print('../../figs/ss_cutbacks_vs_capacities_pressure_ratios','-dpng','-r300');


figure(9)
%mean turning angle vs SS cutback amount
plot(0:10:100, turning_angle_means(1:11)-turning_angle_means(1,1), 'ko-')

xlabel('Cutback amount, %')
ylabel('\Delta turning angle, degrees')

%Here we set up the axes
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) figure_width*100, figure_height*100]); %<- Set size
set(gca, 'FontSize', font_size, 'LineWidth', axes_line_width); %<- Set properties
set(gca,'FontName','Charter','FontSize',font_size)
% Here we preserve the size of the image when we save it.
set(gcf,'InvertHardcopy','on');
set(gcf,'PaperUnits', 'inches');
papersize = get(gcf, 'PaperSize');
left = (papersize(1)- figure_width)/2;
bottom = (papersize(2)- figure_height)/2;
myfiguresize = [left, bottom, figure_width, figure_height];
set(gcf,'PaperPosition', myfiguresize);
print('../../figs/ss_cutbacks_vs_turning_angles','-dpng','-r300');


figure(10)
%mean loss vs SS cutback amount
plot(0:10:100, losses_tp_means(1:11)*100, 'ro-')
hold on
plot(0:10:100, losses_ke_means(1:11)*100, 'bo-')

xlabel('Cutback amount, %')
ylabel('Overall loss at domain outlet, %')
legend('Total pressure loss', 'Kinetic energy loss', 'Location', 'northwest')

%Here we set up the axes
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) figure_width*100, figure_height*100]); %<- Set size
set(gca, 'FontSize', font_size, 'LineWidth', axes_line_width); %<- Set properties
set(gca,'FontName','Charter','FontSize',font_size)
% Here we preserve the size of the image when we save it.
set(gcf,'InvertHardcopy','on');
set(gcf,'PaperUnits', 'inches');
papersize = get(gcf, 'PaperSize');
left = (papersize(1)- figure_width)/2;
bottom = (papersize(2)- figure_height)/2;
myfiguresize = [left, bottom, figure_width, figure_height];
set(gcf,'PaperPosition', myfiguresize);
print('../../figs/ss_cutbacks_vs_losses','-dpng','-r300');

figure(11)
%mean loss vs mean capacity delta
plot(100*T900_TE_capacities_2D_oversampled(3354,2:12)./T900_TE_capacities_2D_oversampled(3354,2)-100, 100*losses_tp_means(1:11)-100*losses_tp_means(1,1), 'ro-')
hold on
plot(100*T900_TE_capacities_2D_oversampled(3354,2:12)./T900_TE_capacities_2D_oversampled(3354,2)-100, 100*losses_ke_means(1:11)-100*losses_ke_means(1,1), 'bo-')

xlabel('\Delta capacity, %')
ylabel('\Delta loss at domain outlet, %')
legend('Total pressure loss', 'Kinetic energy loss', 'Location', 'northwest')

%Here we set up the axes
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) figure_width*100, figure_height*100]); %<- Set size
set(gca, 'FontSize', font_size, 'LineWidth', axes_line_width); %<- Set properties
set(gca,'FontName','Charter','FontSize',font_size)
% Here we preserve the size of the image when we save it.
set(gcf,'InvertHardcopy','on');
set(gcf,'PaperUnits', 'inches');
papersize = get(gcf, 'PaperSize');
left = (papersize(1)- figure_width)/2;
bottom = (papersize(2)- figure_height)/2;
myfiguresize = [left, bottom, figure_width, figure_height];
set(gcf,'PaperPosition', myfiguresize);
print('../../figs/ss_capacities_vs_losses','-dpng','-r300');

figure(12)
%turning angle surveys
plot(tp00(:,1)/(tp00(17,1)-tp00(1,1)),turning_angles_00(:)-turning_angles_00(9), 'k.-')
hold on
plot(tp02(:,1)/(tp00(17,1)-tp00(1,1)),turning_angles_02(:)-turning_angles_00(9), 'r.-')
plot(tp04(:,1)/(tp00(17,1)-tp00(1,1)),turning_angles_04(:)-turning_angles_00(9), 'g.-')
plot(tp06(:,1)/(tp00(17,1)-tp00(1,1)),turning_angles_06(:)-turning_angles_00(9), 'b.-')
plot(tp08(:,1)/(tp00(17,1)-tp00(1,1)),turning_angles_08(:)-turning_angles_00(9), 'c.-')
plot(tp10(:,1)/(tp00(17,1)-tp00(1,1)),turning_angles_10(:)-turning_angles_00(9), 'm.-')

xlabel('Distance along outlet plane')
ylabel('\Delta turning angle, degrees')

%Here we set up the axes
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) figure_width*100, figure_height*100]); %<- Set size
set(gca, 'FontSize', font_size, 'LineWidth', axes_line_width); %<- Set properties
set(gca,'FontName','Charter','FontSize',font_size)
% Here we preserve the size of the image when we save it.
set(gcf,'InvertHardcopy','on');
set(gcf,'PaperUnits', 'inches');
papersize = get(gcf, 'PaperSize');
left = (papersize(1)- figure_width)/2;
bottom = (papersize(2)- figure_height)/2;
myfiguresize = [left, bottom, figure_width, figure_height];
set(gcf,'PaperPosition', myfiguresize);
print('../../figs/ss_cutbacks_turning_angle_surveys','-dpng','-r300');









%%%%%%%%%%%%%%%%%%%%%%%% PLOTS for PS CUTBACKS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



figure(13)
%novel TEs, just the reruns

%plot(T900_capacities_2D(18:28,1), 100*T900_capacities_2D(18:28,2)./T900_capacities_2D(18,2), 'b.')
plot(T900_novel_capacities_2D(:,1), 100*T900_novel_capacities_2D(:,2)./T900_capacities_2D(18,2), 'k.')
hold on
plot(T900_novel_capacities_2D(:,1), 100*T900_novel_capacities_2D(:,10)./T900_capacities_2D(18,2), 'r.')
plot(T900_novel_capacities_2D(:,1), 100*T900_novel_capacities_2D(:,9)./T900_capacities_2D(18,2), 'g.')
plot(T900_novel_capacities_2D(:,1), 100*T900_novel_capacities_2D(:,8)./T900_capacities_2D(18,2), 'b.')
plot(T900_novel_capacities_2D(:,1), 100*T900_novel_capacities_2D(:,6)./T900_capacities_2D(18,2), 'c.')

%plot(T900_capacities_2D_oversampled(:,1), 100*T900_capacities_2D_oversampled(:,2)./T900_capacities_2D(18,2), 'k-')
plot(T900_novel_capacities_2D_oversampled(:,1), 100*T900_novel_capacities_2D_oversampled(:,2)./T900_capacities_2D(18,2), 'k-')
plot(T900_novel_capacities_2D_oversampled(:,1), 100*T900_novel_capacities_2D_oversampled(:,10)./T900_capacities_2D(18,2), 'r-')
plot(T900_novel_capacities_2D_oversampled(:,1), 100*T900_novel_capacities_2D_oversampled(:,9)./T900_capacities_2D(18,2), 'g-')
plot(T900_novel_capacities_2D_oversampled(:,1), 100*T900_novel_capacities_2D_oversampled(:,8)./T900_capacities_2D(18,2), 'b-')
plot(T900_novel_capacities_2D_oversampled(:,1), 100*T900_novel_capacities_2D_oversampled(:,6)./T900_capacities_2D(18,2), 'c-')

xlim([1.5 2.5])
xlabel('NGV pressure ratio (p_{01}/p_{2})')
ylabel('\Delta capacity, %')
design_line = xline(1.79,'-k',{'Design'},'FontName','Charter','FontSize',font_size);
design_line.LabelVerticalAlignment = 'bottom';
legend('0%', '25%', '50%', '75%', '100%', 'Location', 'southeast')

%Here we set up the axes
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) figure_width*100, figure_height*100]); %<- Set size
set(gca, 'FontSize', font_size, 'LineWidth', axes_line_width); %<- Set properties
set(gca,'FontName','Charter','FontSize',font_size)
% Here we preserve the size of the image when we save it.
set(gcf,'InvertHardcopy','on');
set(gcf,'PaperUnits', 'inches');
papersize = get(gcf, 'PaperSize');
left = (papersize(1)- figure_width)/2;
bottom = (papersize(2)- figure_height)/2;
myfiguresize = [left, bottom, figure_width, figure_height];
set(gcf,'PaperPosition', myfiguresize);
print('../../figs/ps_cutbacks_capacities_trends','-dpng','-r300');





figure(14)
plot(25*ps_cutback_throat_widths(:,1), 100*ps_cutback_throat_widths(:,2)/ps_cutback_throat_widths(1,2)-100, 'ko-')
xlabel('Cutback amount, %')
ylabel('Change in throat width, %')
%Here we set up the axes
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) figure_width*100, figure_height*100]); %<- Set size
set(gca, 'FontSize', font_size, 'LineWidth', axes_line_width); %<- Set properties
set(gca,'FontName','Charter','FontSize',font_size)
% Here we preserve the size of the image when we save it.
set(gcf,'InvertHardcopy','on');
set(gcf,'PaperUnits', 'inches');
papersize = get(gcf, 'PaperSize');
left = (papersize(1)- figure_width)/2;
bottom = (papersize(2)- figure_height)/2;
myfiguresize = [left, bottom, figure_width, figure_height];
set(gcf,'PaperPosition', myfiguresize);
print('../../figs/ps_cutbacks_vs_throat_widths','-dpng','-r300');


figure(16)
plot(0:25:100, 100*ps_cutback_capacities_design(:,1)/ps_cutback_capacities_design(1,1)-100, 'ro-')
hold on
plot(0:25:100, 100*ps_cutback_capacities_design(:,2)/ps_cutback_capacities_design(1,2)-100, 'bo-')
plot(0:25:100, 100*ps_cutback_capacities_design(:,3)/ps_cutback_capacities_design(1,3)-100, 'mo-')
xlabel('Cutback amount, %')
ylabel('\Delta Capacity, %')
legend('PR = 1.50', 'PR = 1.79 (design)', 'PR = 3.33', 'Location', 'southeast')
%ylim([95 104])
%Here we set up the axes
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) figure_width*100, figure_height*100]); %<- Set size
set(gca, 'FontSize', font_size, 'LineWidth', axes_line_width); %<- Set properties
set(gca,'FontName','Charter','FontSize',font_size)
% Here we preserve the size of the image when we save it.
set(gcf,'InvertHardcopy','on');
set(gcf,'PaperUnits', 'inches');
papersize = get(gcf, 'PaperSize');
left = (papersize(1)- figure_width)/2;
bottom = (papersize(2)- figure_height)/2;
myfiguresize = [left, bottom, figure_width, figure_height];
set(gcf,'PaperPosition', myfiguresize);
print('../../figs/ps_cutbacks_vs_capacities','-dpng','-r300');


figure(17)
%mean turning angle vs PS cutback amount
plot(0:25:100, turning_angle_ps_means(1:5)-turning_angle_ps_means(1,1), 'ko-')

xlabel('Cutback amount, %')
ylabel('\Delta turning angle, degrees')

%Here we set up the axes
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) figure_width*100, figure_height*100]); %<- Set size
set(gca, 'FontSize', font_size, 'LineWidth', axes_line_width); %<- Set properties
set(gca,'FontName','Charter','FontSize',font_size)
% Here we preserve the size of the image when we save it.
set(gcf,'InvertHardcopy','on');
set(gcf,'PaperUnits', 'inches');
papersize = get(gcf, 'PaperSize');
left = (papersize(1)- figure_width)/2;
bottom = (papersize(2)- figure_height)/2;
myfiguresize = [left, bottom, figure_width, figure_height];
set(gcf,'PaperPosition', myfiguresize);
print('../../figs/ps_cutbacks_vs_turning_angles','-dpng','-r300');


figure(18)
%mean loss vs PS cutback amount
plot(0:25:100, losses_tp_ps_means(1:5)*100, 'ro-')
hold on
plot(0:25:100, losses_ke_ps_means(1:5)*100, 'bo-')

xlabel('Cutback amount, %')
ylabel('Overall loss at domain outlet, %')
legend('Total pressure loss', 'Kinetic energy loss', 'Location', 'best')

%Here we set up the axes
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) figure_width*100, figure_height*100]); %<- Set size
set(gca, 'FontSize', font_size, 'LineWidth', axes_line_width); %<- Set properties
set(gca,'FontName','Charter','FontSize',font_size)
% Here we preserve the size of the image when we save it.
set(gcf,'InvertHardcopy','on');
set(gcf,'PaperUnits', 'inches');
papersize = get(gcf, 'PaperSize');
left = (papersize(1)- figure_width)/2;
bottom = (papersize(2)- figure_height)/2;
myfiguresize = [left, bottom, figure_width, figure_height];
set(gcf,'PaperPosition', myfiguresize);
print('../../figs/ps_cutbacks_vs_losses','-dpng','-r300');






figure(19)
%mean loss vs mean capacity delta
plot(100*T900_ps_capacities_for_plotting(1,:)./T900_ps_capacities_for_plotting(1,1)-100, 100*losses_tp_ps_means(:)-100*losses_tp_ps_means(1,1), 'ro-')
hold on
plot(100*T900_ps_capacities_for_plotting(1,:)./T900_ps_capacities_for_plotting(1,1)-100, 100*losses_ke_ps_means(:)-100*losses_ke_ps_means(1,1), 'bo-')




xlabel('\Delta capacity, %')
ylabel('\Delta loss at domain outlet, %')
legend('Total pressure loss', 'Kinetic energy loss', 'Location', 'northwest')

%Here we set up the axes
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) figure_width*100, figure_height*100]); %<- Set size
set(gca, 'FontSize', font_size, 'LineWidth', axes_line_width); %<- Set properties
set(gca,'FontName','Charter','FontSize',font_size)
% Here we preserve the size of the image when we save it.
set(gcf,'InvertHardcopy','on');
set(gcf,'PaperUnits', 'inches');
papersize = get(gcf, 'PaperSize');
left = (papersize(1)- figure_width)/2;
bottom = (papersize(2)- figure_height)/2;
myfiguresize = [left, bottom, figure_width, figure_height];
set(gcf,'PaperPosition', myfiguresize);
print('../../figs/ps_capacities_vs_losses','-dpng','-r300');


figure(20)
%PS turning angle surveys
plot(tpps0(:,1)/(tpps0(17,1)-tpps0(1,1)),turning_angles_ps_0(:)-turning_angles_ps_0(9), 'k.-')
hold on
plot(tpps1(:,1)/(tpps1(17,1)-tpps1(1,1)),turning_angles_ps_1(:)-turning_angles_ps_0(9), 'r.-')
plot(tpps2(:,1)/(tpps2(17,1)-tpps2(1,1)),turning_angles_ps_2(:)-turning_angles_ps_0(9), 'g.-')
plot(tpps3(:,1)/(tpps3(17,1)-tpps3(1,1)),turning_angles_ps_3(:)-turning_angles_ps_0(9), 'b.-')
plot(tpps4(:,1)/(tpps4(17,1)-tpps4(1,1)),turning_angles_ps_4(:)-turning_angles_ps_0(9), 'c.-')

xlabel('Distance along outlet plane')
ylabel('\Delta turning angle, degrees')

%Here we set up the axes
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) figure_width*100, figure_height*100]); %<- Set size
set(gca, 'FontSize', font_size, 'LineWidth', axes_line_width); %<- Set properties
set(gca,'FontName','Charter','FontSize',font_size)
% Here we preserve the size of the image when we save it.
set(gcf,'InvertHardcopy','on');
set(gcf,'PaperUnits', 'inches');
papersize = get(gcf, 'PaperSize');
left = (papersize(1)- figure_width)/2;
bottom = (papersize(2)- figure_height)/2;
myfiguresize = [left, bottom, figure_width, figure_height];
set(gcf,'PaperPosition', myfiguresize);
print('../../figs/ps_cutbacks_turning_angle_surveys','-dpng','-r300');
















