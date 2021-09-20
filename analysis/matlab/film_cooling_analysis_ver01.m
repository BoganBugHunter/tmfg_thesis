clear variables;
close all;
load workspaces/sch_data.mat;


%finding the plenum total pressure from a given mass flow rate fraction

%physical properties
R = 287.058; %[Joules/(kg.K)]
cp = 1121; %[Joules/(kg.K)]
cv = 834; %[Joules/(kg.K)]
gamma = cp/cv;

fluent_operating_pressure = 4285900;

p_m_benchmark = -1e6 + fluent_operating_pressure;

coolant_ratio = 0.01;

%2D hole width
A_h = 7.4993e-4;

%inlet total pressure
p_0m = 4325590;

%total temperatures
T_0m = 300;
T_0c = 300;

%correction between 2D and Q3D
q3d_strip_width = 0.0029;
ratio_of_slot_to_hole = (4/pi)*q3d_strip_width/A_h;



%use the benchmark coolant MFR to work out the plenum total pressure,
%iteratively

m_m_benchmark = 121;
m_c_benchmark = m_m_benchmark*coolant_ratio;
zeta_p = 0.8;

number_of_iterations = 10;

%array to store and plot the convergence for both possible methods
array = zeros(number_of_iterations,1);

%initial guess
p_0c_benchmark = p_m_benchmark;

for i = 1:number_of_iterations
  p_0c_benchmark = p_m_benchmark*( 1 - ( R*(gamma-1)*T_0c*m_c_benchmark^2*p_0c_benchmark^((2-2*gamma)/gamma) )/( 2*gamma*A_h^2*p_m_benchmark^(2/gamma) ) )^(gamma/(1-gamma));
  array(i,1) = p_0c_benchmark;
end

ratio_of_coolant_to_main_tp = p_0c_benchmark/p_0m;

resulting_coolant_mass_flow_rate = (p_0c_benchmark/sqrt(T_0c))*A_h*sqrt(gamma/R)*(p_m_benchmark/p_0c_benchmark)^(1/gamma)*sqrt( (2/(gamma-1))*( 1 - (p_m_benchmark/p_0c_benchmark)^((gamma-1)/gamma) ) );
resulting_coolant_ratio = resulting_coolant_mass_flow_rate/m_m_benchmark;



%to find the correct plenum total pressure
p_0cp = ( p_0c_benchmark - zeta_p*p_m_benchmark )/( 1 - zeta_p );
ratio_of_plenum_to_main_tp = p_0cp/p_0m;

p_0cp_gauge = p_0cp - fluent_operating_pressure;







%plot the CFD data%%%%%%%



%mass_flow_rates columns:
%1) hole position (not normalised)
%2) inlet MFR
%3) outlet MFR
%4) plenum MFR
%5) local surface static pressure (gauge in Fluent)



%outlet_total_pressures_averages columns:
%1) hole position (not normalised)
%2) area-averaged outlet total pressure (gauge in Fluent)



%normalise the hole positions

hole_positions = (50/27)*(mass_flow_rates(:,1) - 26);


%work out the theoretical hole exit total pressures and thus capacities

m_c = mass_flow_rates(:,4);

p_m = mass_flow_rates(:,5) + fluent_operating_pressure;

p_0c = zeta_p*p_m(:) + (1-zeta_p)*p_0cp;

G_c = ( sqrt(T_0c)./p_0c )*m_c_benchmark;


%work out the inlet capacities

m_in = mass_flow_rates(:,2);

G_in = ( sqrt(T_0m)./p_0m ).*m_in;

G_out = G_c + G_in;



%define a surface pressure coefficient
C_p = (p_m-min(p_m))/(max(p_m)-min(p_m));

%define a surface isentropic Mach number
M_surf = sqrt( (2/(gamma-1))*( (p_0m./p_m).^((gamma-1)/gamma) -1 ) );

p_2 = 2625200; %outlet static pressure

%define a loss for each hole position
total_pressures = outlet_total_pressure_averages(:,2);
total_pressures_absolute = total_pressures + fluent_operating_pressure;
loss = (p_0m - total_pressures_absolute)/(p_0m - p_2);





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


geometric_throat_location = 38;
geometric_throat_normalised = (50/27)*(geometric_throat_location - 26);


figure(1)
plot(hole_positions,100*(G_in/G_in(1,1) - 1), 'r-o')
hold on
plot(hole_positions,100*(G_out/G_out(1,1) - 1), 'b-o')

xlabel('Distance along suction side, % of range')
ylabel('\Delta capacity, %')
throat_line = xline(geometric_throat_normalised,'-k',{'Geometric throat'},'FontName','Charter','FontSize',font_size);
throat_line.LabelVerticalAlignment = 'top';
legend('\Gamma_{in}', '\Gamma_{out}', 'Location', 'southeast')

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
print('../../figs/sch_hole_location_vs_capacity','-dpng','-r300');


figure(2)
plot(hole_positions,100*(G_c/G_c(1,1) - 1), 'k-o')

xlabel('Distance along suction side, % of range')
ylabel('\Delta  coolant capacity, %')
throat_line = xline(geometric_throat_normalised,'-k',{'Geometric throat'},'FontName','Charter','FontSize',font_size);
throat_line.LabelVerticalAlignment = 'bottom';

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
print('../../figs/sch_hole_location_vs_coolant_capacity','-dpng','-r300');


figure(3)
plot(  C_p,  100*(G_in/G_in(1,1) - 1), 'r-o')
hold on
plot(  C_p,  100*(G_out/G_out(1,1) - 1), 'b-o')

xlabel('Surface static pressure coefficient')
ylabel('\Delta  capacity, %')

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
print('../../figs/sch_static_pressure_coefficient_vs_capacity','-dpng','-r300');


figure(4)
plot(  M_surf,  100*(G_in/G_in(1,1) - 1), 'r-o')
hold on
plot(  M_surf,  100*(G_out/G_out(1,1) - 1), 'b-o')

xlabel('Surface isentropic Mach number')
ylabel('\Delta capacity, %')

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
print('../../figs/sch_surface_isentropic_mach_number_vs_capacity','-dpng','-r300');


figure(5)
plot(  C_p,  100*(G_c/G_c(1,1) - 1), 'k-o')

xlabel('Surface static pressure coefficient')
ylabel('\Delta  coolant capacity, %')

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
print('../../figs/sch_static_pressure_coefficient_vs_coolant_capacity','-dpng','-r300');


figure(6)
plot(hole_positions,100*loss, 'k-o')

xlabel('Distance along suction side, % of range')
ylabel('Overall loss at domain outlet, %')
throat_line = xline(geometric_throat_normalised,'-k',{'Geometric throat'},'FontName','Charter','FontSize',font_size);
throat_line.LabelVerticalAlignment = 'bottom';

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
print('../../figs/sch_hole_location_vs_loss','-dpng','-r300');












