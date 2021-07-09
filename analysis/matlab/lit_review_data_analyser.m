close all
clear variables
load lit_review_data.mat


%physical properties
R = 287.058; %[Joules/(kg.K)]
cp = 1121; %[Joules/(kg.K)]
cv = 834; %[Joules/(kg.K)]
gamma = cp/cv; 


figure(1)
hold on
%axis equal
xlim([0 3])
ylim([0 8])
plot(raffel122coolant0(:,1),raffel122coolant0(:,2));
plot(raffel122measuredtotal(:,1),raffel122measuredtotal(:,2));
plot(raffel122upstreamofte(:,1),raffel122upstreamofte(:,2));
plot(raffel122teloss(:,1),raffel122teloss(:,2), 'k--');

figure(2)
hold on
%axis equal
xlim([0 3])
ylim([0 8])
plot(raffel134coolant0(:,1),raffel134coolant0(:,2));
plot(raffel134measuredtotal(:,1),raffel134measuredtotal(:,2));
plot(raffel134upstreamofte(:,1),raffel134upstreamofte(:,2));
plot(raffel134teloss(:,1),raffel134teloss(:,2), 'k--');

figure(3)
hold on
%axis equal
xlim([0.1 0.8])
ylim([0.02 0.2])
plot(backdatrinidade_loss_models(:,1),backdatrinidade_loss_models(:,2), 'k*')
plot(backdatrinidade_loss_models(:,3),backdatrinidade_loss_models(:,4), 'ko')
plot(backdatrinidade_loss_models(:,5),backdatrinidade_loss_models(:,6), 'k.')

figure(4)
hold on
%axis equal
xlim([-0.005 0.06])
ylim([0.012 0.028])
plot(gaoTE_cpt_and_alpha_vs_blowing(:,1),gaoTE_cpt_and_alpha_vs_blowing(:,2), 'k*-')
plot(gaoTE_cpt_and_alpha_vs_blowing(:,3),gaoTE_cpt_and_alpha_vs_blowing(:,4), 'ko-')



%compare Giel's and Gao's loss definitions

PR_critical = ((gamma+1)/2)^(gamma/(1-gamma));

p01 = 1;
p2 = p01*PR_critical;
p0c = 0.2*p01;

Delta_TP = 0:0.001:0.4*(p01-p2);
p02 = p01 - Delta_TP + p0c;

CpGiel = (p01-p02)./(p01-p2);
CpGao = (p01-p02)./(p02-p2);

CpGaoSeries1 = CpGiel;
CpGaoSeries2 = CpGiel + CpGiel.^2;
CpGaoSeries3 = CpGiel + CpGiel.^2 + CpGiel.^3;
CpGaoSeries4 = CpGiel + CpGiel.^2 + CpGiel.^3 + CpGiel.^4;
CpGaoSeries5 = CpGiel + CpGiel.^2 + CpGiel.^3 + CpGiel.^4 + CpGiel.^5;

CpGielSeries1 = CpGao;
CpGielSeries2 = CpGao - CpGao.^2;
CpGielSeries3 = CpGao - CpGao.^2 + CpGao.^3;
CpGielSeries4 = CpGao - CpGao.^2 + CpGao.^3 - CpGao.^4;
CpGielSeries5 = CpGao - CpGao.^2 + CpGao.^3 - CpGao.^4 + CpGao.^5;


figure(5)
plot(Delta_TP,CpGaoSeries1, 'm-');
hold on
plot(Delta_TP,CpGaoSeries2, 'c-');
plot(Delta_TP,CpGaoSeries3, 'r-');
plot(Delta_TP,CpGaoSeries4, 'b-');
plot(Delta_TP,CpGaoSeries5, 'g-');
plot(Delta_TP,CpGao, 'k--');

figure(6)
plot(Delta_TP,CpGielSeries1, 'm-');
hold on
plot(Delta_TP,CpGielSeries2, 'c-');
plot(Delta_TP,CpGielSeries3, 'r-');
plot(Delta_TP,CpGielSeries4, 'b-');
plot(Delta_TP,CpGielSeries5, 'g-');
plot(Delta_TP,CpGiel, 'k--');

