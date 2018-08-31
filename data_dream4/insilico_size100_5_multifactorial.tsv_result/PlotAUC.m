addpath('/Users/chenx/Documents/Exp/bLARS-Exp/ARNI');
score_mat = csvread('idx_aupr.csv');
dataname = 'DREAM4 Multifactorial Network5';
texname=  'DREAM4_Multifactorial_Network_5_AUPR';

changeplot;
t = 1:1:5;
y1 = score_mat(1,t);
y2 = score_mat(2,t);
y3 = score_mat(3,t);
y4 = score_mat(4,t);
y5 = score_mat(5,t);

grid on;
hold on;
plot(t,y1);
plot(t,y2);
plot(t,y3);
plot(t,y4);
plot(t,y5);
hold off;

legend('BS=100','BS=200','BS=300','BS=400','BS=500')
title(dataname);
xlabel('Number of ARNI steps');
ylabel('AUPR');

saveas(gcf,[pwd '/aupr_vs_S.png']);
matlab2tikz([texname '.tikz'],'height', '\fheight', 'width', '\fwidth');

close all;
clear all;

score_mat = csvread('idx_auroc.csv');
dataname = 'DREAM4 Multifactorial Network5';
texname=  'DREAM4_Multifactorial_Network_5_AUROC';

changeplot;
t = 1:1:5;
y1 = score_mat(1,t);
y2 = score_mat(2,t);
y3 = score_mat(3,t);
y4 = score_mat(4,t);
y5 = score_mat(5,t);

grid on;
hold on;
plot(t,y1);
plot(t,y2);
plot(t,y3);
plot(t,y4);
plot(t,y5);
hold off;

legend('BS=100','BS=200','BS=300','BS=400','BS=500')
title(dataname);
xlabel('Number of ARNI steps');
ylabel('AUROC');

saveas(gcf,[pwd '/auroc_vs_S.png']);
matlab2tikz([texname '.tikz'],'height', '\fheight', 'width', '\fwidth');

close all;