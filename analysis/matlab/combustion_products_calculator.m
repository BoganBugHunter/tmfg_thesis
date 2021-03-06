close all;
clear variables;


%physical constants for air at stp
cp_air = 1010; %[Joules/(kg.K)]
cv_air = 718; %[Joules/(kg.K)]
R_air = cp_air-cv_air; %[Joules/(kg.K)]
g_air = cp_air/cv_air; 

%physical constants for nitrogen
cp_n2 = 1040; %[Joules/(kg.K)]
cv_n2 = 743; %[Joules/(kg.K)]
R_n2 = cp_n2-cv_n2; %[Joules/(kg.K)]
g_n2 = cp_n2/cv_n2;
m_n2 = 28; %[g/mol]

%physical constants for oxygen
cp_o2 = 919; %[Joules/(kg.K)]
cv_o2 = 659; %[Joules/(kg.K)]
R_o2 = cp_o2-cv_o2; %[Joules/(kg.K)]
g_o2 = cp_o2/cv_o2;
m_o2 = 32; %[g/mol]

%physical constants for C02
cp_co2 = 844; %[Joules/(kg.K)]
cv_co2 = 655; %[Joules/(kg.K)]
R_co2 = cp_co2-cv_co2; %[Joules/(kg.K)]
g_co2 = cp_co2/cv_co2;
m_co2 = 44; %[g/mol]

%physical constants for steam at as close to STP while still superheated
cp_h2o = 2175; %[Joules/(kg.K)]
cv_h2o = 1616; %[Joules/(kg.K)]
R_h2o = cp_h2o-cv_h2o; %[Joules/(kg.K)]
g_h2o = cp_h2o/cv_h2o;
m_h2o = 18; %[g/mol]

%properties of air as calculated from its stated constituents
cp_air_calc = (76.70*cp_n2 + 23.30*cp_o2)/100; %[Joules/(kg.K)]
cv_air_calc = (76.70*cv_n2 + 23.30*cv_o2)/100; %[Joules/(kg.K)]
R_air_calc = cp_air_calc-cv_air_calc; %[Joules/(kg.K)]
g_air_calc = cp_air_calc/cv_air_calc;

%properties of combustion product as calculated from its stated
%constituents
cp_air_prod = (71.82*cp_n2 + 20.00*cp_co2 + 8.18*cp_h2o)/100; %[Joules/(kg.K)]
cv_air_prod = (71.82*cv_n2 + 20.00*cv_co2 + 8.18*cv_h2o)/100; %[Joules/(kg.K)]
R_air_prod = cp_air_prod-cv_air_prod; %[Joules/(kg.K)]
g_air_prod = cp_air_prod/cv_air_prod;

%properties of hot air (2000K)
cp_air_hot = 1241; %[Joules/(kg.K)]
cv_air_hot = 953.5; %[Joules/(kg.K)]
R_air_hot = cp_air_hot-cv_air_hot; %[Joules/(kg.K)]
g_air_hot = cp_air_hot/cv_air_hot;

%how far off are we with the calculated air?
percentage_change_in_cp_for_calc = 100*(cp_air_calc/cp_air - 1);
percentage_change_in_cv_for_calc = 100*(cv_air_calc/cv_air - 1);
percentage_change_in_R_for_calc = 100*(R_air_calc/R_air - 1);
percentage_change_in_g_for_calc = 100*(g_air_calc/g_air - 1);

%how far off are we with the combustion products?
percentage_change_in_cp_for_prod = 100*(cp_air_prod/cp_air - 1);
percentage_change_in_cv_for_prod = 100*(cv_air_prod/cv_air - 1);
percentage_change_in_R_for_prod = 100*(R_air_prod/R_air - 1);
percentage_change_in_g_for_prod = 100*(g_air_prod/g_air - 1);

%how far off are we with the hot air?
percentage_change_in_cp_for_hot = 100*(cp_air_hot/cp_air - 1);
percentage_change_in_cv_for_hot = 100*(cv_air_hot/cv_air - 1);
percentage_change_in_R_for_hot = 100*(R_air_hot/R_air - 1);
percentage_change_in_g_for_hot = 100*(g_air_hot/g_air - 1);


%work out the 1D capacity trends for STP air, combustion products, and hot
%air

%the range of pressure ratios
r = 0.001:0.001:1;

%the capacity trends

Gamma_air = sqrt( 2*g_air/(R_air*(g_air-1))*r.^(2/g_air).*(1-r.^((g_air-1)...
    /g_air)) );
r_crit_air = ((g_air+1)/2)^(g_air/(1-g_air));
Gamma_crit_air = sqrt(g_air/R_air)*((g_air+1)/2)^((1+g_air)/(2*(1-g_air)));
for i = 1:length(r)
    if r(i) < r_crit_air
        Gamma_air(i) = Gamma_crit_air;
    end
end

Gamma_air_calc = sqrt( 2*g_air_calc/(R_air_calc*(g_air_calc-1))*r.^(2/...
    g_air_calc).*(1-r.^((g_air_calc-1)/g_air_calc)) );
r_crit_air_calc = ((g_air_calc+1)/2)^(g_air_calc/(1-g_air_calc));
Gamma_crit_air_calc = sqrt(g_air_calc/R_air_calc)*((g_air_calc+1)/2)^...
    ((1+g_air_calc)/(2*(1-g_air_calc)));
for i = 1:length(r)
    if r(i) < r_crit_air_calc
        Gamma_air_calc(i) = Gamma_crit_air_calc;
    end
end

Gamma_air_prod = sqrt( 2*g_air_prod/(R_air_prod*(g_air_prod-1))*r.^(2/...
    g_air_prod).*(1-r.^((g_air_prod-1)/g_air_prod)) );
r_crit_air_prod = ((g_air_prod+1)/2)^(g_air_prod/(1-g_air_prod));
Gamma_crit_air_prod = sqrt(g_air_prod/R_air_prod)*((g_air_prod+1)/2)^...
    ((1+g_air_prod)/(2*(1-g_air_prod)));
for i = 1:length(r)
    if r(i) < r_crit_air_prod
        Gamma_air_prod(i) = Gamma_crit_air_prod;
    end
end

Gamma_air_hot = sqrt( 2*g_air_hot/(R_air_hot*(g_air_hot-1))*r.^...
    (2/g_air_hot).*(1-r.^((g_air_hot-1)/g_air_hot)) );
r_crit_air_hot = ((g_air_hot+1)/2)^(g_air_hot/(1-g_air_hot));
Gamma_crit_air_hot = sqrt(g_air_hot/R_air_hot)*((g_air_hot+1)/2)^...
    ((1+g_air_hot)/(2*(1-g_air_hot)));
for i = 1:length(r)
    if r(i) < r_crit_air_hot
        Gamma_air_hot(i) = Gamma_crit_air_hot;
    end
end

%normalise each trend against the cold air trend
Factor_air_calc = 100*(Gamma_air_calc./Gamma_air-1);
Factor_air_prod = 100*(Gamma_air_prod./Gamma_air-1);
Factor_air_hot = 100*(Gamma_air_hot./Gamma_air-1);

%work out the inverse pressure ratios and plot the capacities the other way
%round
r_inv = 1./r;


% Default figure settings
figure_width = 4;     % Width in inches
figure_height = 4;    % Height in inches
axes_line_width = 0.75;    % AxesLineWidth
font_size = 14;      % Fontsize
line_width = 1.5;      % LineWidth
marker_size = 8;       % MarkerSize

figure(1)
plot(r_inv,100*(Gamma_air/Gamma_crit_air-1),'k')

hold on
plot(r_inv,100*(Gamma_air_hot/Gamma_crit_air-1),'r')
plot(r_inv,100*(Gamma_air_prod/Gamma_crit_air-1),'b')
xline(1/r_crit_air,'k')
xline(1/r_crit_air_hot,'r')
xline(1/r_crit_air_prod,'b')
xlim([1.5,2])
legend('Air, 300 K', 'Air, 2000 K', 'Product', 'Location', 'Southeast')
xlabel('NGV pressure ratio, $\frac{p_{01}}{p_2}$','Interpreter','latex')
ylabel('\Delta capacity, %')
%Here we set up the axes
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) figure_width*100, figure_height*100]);
set(gca, 'FontSize', font_size, 'LineWidth', axes_line_width);
set(gca,'FontName','Charter','FontSize',font_size)
% Here we preserve the size of the image when we save it.
set(gcf,'InvertHardcopy','on');
set(gcf,'PaperUnits', 'inches');
papersize = get(gcf, 'PaperSize');
left = (papersize(1)- figure_width)/2;
bottom = (papersize(2)- figure_height)/2;
myfiguresize = [left, bottom, figure_width, figure_height];
set(gcf,'PaperPosition', myfiguresize);
print('../../figs/combustion_products_capacities','-dpng','-r300');

figure(2)
plot(r_inv,Factor_air_prod,'b')
hold on
plot(r_inv,Factor_air_hot,'r')
xlim([1.5,2])
xlabel('NGV pressure ratio, $\frac{p_{01}}{p_2}$','Interpreter','latex')
ylabel('Error in capacity, %')
design_line = xline(1.79,'-k',{'Design'},'FontName','Charter','FontSize',...
    font_size);
design_line.LabelVerticalAlignment = 'top';
design_line.LabelHorizontalAlignment = 'right';
choked_line = xline(1/r_crit_air,'-k',{'Critical'},'FontName','Charter',...
    'FontSize',font_size);
choked_line.LabelVerticalAlignment = 'top';
choked_line.LabelHorizontalAlignment = 'right';
legend('Air, 2000 K', 'Product', 'Location', 'Southwest')
%Here we set up the axes
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) figure_width*100, figure_height*100]);
set(gca, 'FontSize', font_size, 'LineWidth', axes_line_width);
set(gca,'FontName','Charter','FontSize',font_size)
% Here we preserve the size of the image when we save it.
set(gcf,'InvertHardcopy','on');
set(gcf,'PaperUnits', 'inches');
papersize = get(gcf, 'PaperSize');
left = (papersize(1)- figure_width)/2;
bottom = (papersize(2)- figure_height)/2;
myfiguresize = [left, bottom, figure_width, figure_height];
set(gcf,'PaperPosition', myfiguresize);
print('../../figs/combustion_products_capacities_errors','-dpng','-r300');













