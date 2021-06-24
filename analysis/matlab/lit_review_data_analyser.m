close all
clear variables
load lit_review_data.mat

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