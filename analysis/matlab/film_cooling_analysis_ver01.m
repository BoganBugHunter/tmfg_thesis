clear variables;
close all;


%finding the plenum total pressure from a given mass flow rate fraction

%physical properties
R = 287.058; %[Joules/(kg.K)]
cp = 1121; %[Joules/(kg.K)]
cv = 834; %[Joules/(kg.K)]
gamma = cp/cv;

p_w = 3085900;

coolant_ratio = 0.008;

T_0C = 300;

m_dot_main = 120;

m_dot_coolant = m_dot_main*coolant_ratio;

plenum_area = 7.4993e-4;

zeta_p = 0.7;

A = 2*(gamma-1)/gamma;
B = p_w^((gamma-1)/gamma);
C = (1-gamma)/gamma;
D = m_dot_coolant^2/( plenum_area^2*p_w^(2/gamma)*( 2*gamma/(R*T_0C*(gamma-1)) ) );

%iteratively solve for the plenum total pressure

size = 10;

%array to store and plot the convergence for both possible methods
array = zeros(size,1);

%initial guess
p_0c = p_w;

for i = 1:size
   p_0c = ( (1/B)*( 1 - D/(p_0c^A) ) )^(1/C);
   array(i,1) = p_0c;
end

ratio_of_coolant_to_main_tp = p_0c/p_w;



%to find the correct plenum total pressure
p_0cp = ( p_0c - zeta_p*p_w )/( 1 - zeta_p );


ratio_of_plenum_to_main_tp = p_0cp/p_w;



figure(1)
plot(array(:,1))


