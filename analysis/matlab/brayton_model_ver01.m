%code to model the Brayton cycle in a Trent XWB engine for illustrations
% http://dgleich.github.io/hq-matlab-figs/ for good figures

close all;
clear variables;


%physical properties
R = 287.058; %[Joules/(kg.K)]
cp = 1121; %[Joules/(kg.K)]
cv = 834; %[Joules/(kg.K)]
gamma = cp/cv; 
rho_atm = 0.3*1.225; %[kg/m3] at 35000 ft
p_atm = 23800; %[Pascals] at 35000 ft
M_cruise = 0.85;

%the engine's mechanical design
compression_ratio = 1/4; % design constraint
T3 = 2000 + 273.15; %[K] design constraint
m_dot = 1440; %[kg/s] source: newatlas.com/rolls-royce-trent-xwb-engine-a350-xwb/32133

%emergent intitial conditions
p_dyn = M_cruise^2*0.5*gamma*p_atm;
p1 = p_atm + p_dyn;
v1 = 1/rho_atm;
T1 = v1*p1/R;
T1_in_celcius = T1 - 273.15;
s1 = 0; % since the entropy is relative


% 1 to 2 - isentropic compression
v2 = v1*compression_ratio;
p2 = p1*(v2/v1)^(-gamma);
T2 = v2.*p2/R;
s2 = s1;
% and the lines of points
line_12v = linspace(v1,v2,500);
line_12p = p1*(line_12v/v1).^(-gamma);
line_12T = line_12v.*line_12p/R;
line_12s = linspace(s1,s2,500);

% 2 to 3 - heat addition at constant pressure
p3 = p2;
Qin = m_dot*cp*(T3-T2);
Qin_in_hp = Qin*0.001341022;
v3 = R*T3/p3;
s3 = cp*log(T3/T2);
% and the lines of points
line_23p = linspace(p2,p3,500);
line_23v = linspace(v2,v3,500);
line_23T = line_23v.*line_23p/R;
line_23s = cp.*log(line_23T./T2);

% 3 to 4 - isentropic expansion
p4 = p1;
T4 = T3/(p4/p3)^((1-gamma)/gamma);
v4 = R*T4/p4;
s4 = s3;
%and the lines of points
line_34v = linspace(v3,v4,500);
line_34p = p4*(line_34v/v4).^(-gamma);
line_34T = line_34v.*line_34p/R;
line_34s = linspace(s3,s4,500);

% 4 to 1 - heat rejection at constant pressure
line_41v = linspace(v4,v1,500);
line_41p = linspace(p1,p1,500);
line_41T = line_41v.*line_41p/R;
line_41s = cp.*log(line_41T./T1);

% what does the Schlicting skin friction correlation look like as a
% function of Reynolds number?
Re = 0:10000:100e6;
Cf = (2*log10(Re) - 0.65).^-2.3;
Cfroot = 0.5*sqrt(Cf);

%did I calculate y-plus correctly?
mu = 19e-6;
L = 0.08;
p0 = 4.33e6;
T0 = 300;
y = 3e-7;
Re_value = (L/mu)*sqrt(gamma/R)*(p0/sqrt(T0))*(1 + 0.5*(gamma-1))^((gamma+1)/(2*(1-gamma)));
delta = 0.37*L*Re_value^-0.2;
Cf_value = (2*log10(Re_value) - 0.65).^-2.3;
y_plus = (y/mu)*sqrt(gamma/R)*(p0/sqrt(T0))*sqrt(Cf_value/2)*(1 + 0.5*(gamma-1))^((gamma+1)/(2*(1-gamma)));



% Default figure settings
figure_width = 4;     % Width in inches
figure_height = 4;    % Height in inches
axes_line_width = 0.75;    % AxesLineWidth
font_size = 15;      % Fontsize
line_width = 1.5;      % LineWidth
marker_size = 8;       % MarkerSize

%to plot a single p-v diagram when desired
figure (2)
plot(line_12v,line_12p,'k','LineWidth',2);
hold on; % needs to come after the first plot command, otherwise the bounding box is missing
plot(line_23v, line_23p, 'r','LineWidth',2);
plot(line_34v,line_34p,'k','LineWidth',2);
plot(line_41v,line_41p, 'c','LineWidth',2);
text(v1,p1,'1','FontName','Charter','FontSize',font_size,'VerticalAlignment','top','HorizontalAlignment','right')
text(v2,p2,'2','FontName','Charter','FontSize',font_size,'VerticalAlignment','bottom','HorizontalAlignment','right')
text(v3,p3,'3','FontName','Charter','FontSize',font_size,'VerticalAlignment','bottom','HorizontalAlignment','left')
text(v4,p4,'4','FontName','Charter','FontSize',font_size,'VerticalAlignment','top','HorizontalAlignment','left')
xlabel('v');
ylabel('p');
xlim([v2-0.15*(v4-v2) v4+0.15*(v4-v2)]);
ylim([p1-0.15*(p2-p1) p2+0.15*(p2-p1)]);
set(gca,'XTick',[], 'YTick', []) % no numbers for this illustration
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
print('../../figs/brayton_cycle_pv_plot','-dpng','-r300');

%to plot a single T-s diagram when desired
figure (3)
plot(line_12s,line_12T,'k','LineWidth',2);
hold on; % needs to come after the first plot command, otherwise the bounding box is missing
plot(line_23s, line_23T, 'r','LineWidth',2);
plot(line_34s,line_34T,'k','LineWidth',2);
plot(line_41s,line_41T, 'c','LineWidth',2);
text(s1,T1,'1','FontName','Charter','FontSize',font_size,'VerticalAlignment','top','HorizontalAlignment','right')
text(s2,T2,'2','FontName','Charter','FontSize',font_size,'VerticalAlignment','bottom','HorizontalAlignment','right')
text(s3,T3,'3','FontName','Charter','FontSize',font_size,'VerticalAlignment','bottom','HorizontalAlignment','left')
text(s4,T4,'4','FontName','Charter','FontSize',font_size,'VerticalAlignment','top','HorizontalAlignment','left')
xlabel('s');
ylabel('T');
xlim([s1-0.15*(s4-s1) s4+0.15*(s4-s1)]);
ylim([T1-0.15*(T3-T1) T3+0.15*(T3-T1)]);
set(gca,'XTick',[], 'YTick', []) % no numbers for this illustration
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
print('../../figs/brayton_cycle_ts_plot','-dpng','-r300');