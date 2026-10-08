
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% figure1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;clc;
% Locate bundled inputs relative to this saved script, independent of pwd.
figureDir = fileparts(mfilename('fullpath'));
codeRoot = fileparts(figureDir);
dataRoot = fullfile(codeRoot,'figure_data');
addpath(fullfile(codeRoot,'utilities'));
f_file = [fullfile(codeRoot,'figure_outputs','figure1') filesep];
if ~isfolder(f_file)
    mkdir(f_file);
end
ot_color = [226 106 83;242 199 178; 88 118 227;171 198 251]./255;
load(fullfile(dataRoot,'behavioral_results','all_accauc.mat'),'osub','psub');
load(fullfile(dataRoot,'behavioral_results','all_GLM_beta_withGroup_fixExpID.mat'),'beta_ot0','beta_pl0');
load(fullfile(dataRoot,'behavioral_results','sigmoid_data.mat'),'sig_isc','sig_isc0','sig_conf','sig_conf0');

%% figure1D
% synergy(2HBT1)
osc=max([osub.acc_sc(:,1),osub.acc_oc(:,1)],[],2);psc=max([psub.acc_sc(:,1),psub.acc_oc(:,1)],[],2);
onc=max([osub.acc_nc(:,1),osub.acc_ic(:,1)],[],2);pnc=max([psub.acc_nc(:,1),psub.acc_ic(:,1)],[],2);
xx1=[osub.acc_sc(:,2)-osc,psub.acc_sc(:,2)-psc,osub.acc_nc(:,2)-onc,psub.acc_nc(:,2)-pnc];
ax1=figureX1([23 20]);
x_pos=figureTmp1_OT(xx1,3);hold on;plot([xlim],[0 0],'k--','LineWidth',2);
box off;ylabel('synergy effect');set(gca,'XTick',x_pos,'XTickLabel',{'OT','PL','OT','PL'}); gcaf1(gca,35);cutBackColor(gca,gcf);
setFixedPlotArea(ax1, [15 15]);setAxesForPPT(ax1, 35, 2.5);set(gcf,'Color',[1 1 1]);cutBackColor(gca,gcf);
print(gcf,[f_file,'acc_2HBT1.svg'], '-dsvg','-vector','-r600');

%% figure1E 
% chocie sigmoid
ax1=figureX1([23 20]);
hold on
x_tmp=-4:0.01:4;
for i=1:length(sig_isc)
    yv1=predict(sig_isc{i,1},x_tmp');yv2=predict(sig_isc{i,2},x_tmp');
    pl=plot(x_tmp',yv1,"Color",ot_color(1,:));pl.Color(4)=0.1;
    pl=plot(x_tmp',yv2,"Color",ot_color(2,:),'LineStyle','--');pl.Color(4)=0.3;

    yv1=predict(sig_isc{i,3},x_tmp');yv2=predict(sig_isc{i,4},x_tmp');
    pl=plot(x_tmp',yv1,"Color",ot_color(3,:));pl.Color(4)=0.1;
    pl=plot(x_tmp',yv2,"Color",ot_color(4,:),'LineStyle','--');pl.Color(4)=0.1;
end
for i=1:4
plot(x_tmp',predict(sig_isc0{i},x_tmp'),"Color",ot_color(i,:),'LineWidth',3);
end
setFixedPlotArea(ax1, [15 15]);set(gca, 'LineWidth', 3); 
hold off;xlim([-4 4]);ylabel('p (stay)');xlabel('self-other');gcaf1(gca,35);setAxesForPPT(ax1, 35, 2.5);cutBackColor(gca,gcf);
print(gcf,[f_file,'choice_sigmoid.svg'], '-dsvg','-vector','-r600');
% confidence sigmoid
tranlogodds=@(x) 1./(1+exp(-x));
ax1=figureX1([23 20]);
x_tmp=[0:0.01:4]';
hold on
for i=1:length(sig_conf)
    for j=1:4
    yv1=tranlogodds(predict(sig_conf{i,j},x_tmp));
    pl=plot(x_tmp,yv1,"Color",ot_color(j,:));pl.Color(4)=0.1;
    end
end
for j=[3 4 1 2]
plot(x_tmp,tranlogodds(predict(sig_conf0{j},x_tmp)),"Color",ot_color(j,:),'LineWidth',3);
end
setFixedPlotArea(ax1, [15 15]);set(gca, 'LineWidth', 3); 
hold off;xlim([0 4]);ylabel('confidence (2nd)');xlabel({'chosen-unchosen'});gcaf1(gca,35);setAxesForPPT(ax1, 35, 2.5);%cutBackColor(gca,gcf);
print(gcf,[f_file,'confidence_sigmoid0.svg'], '-dsvg','-vector','-r600');


%% choice
% choice ΔV
xx1=[beta_ot0.isc_so1(:,3),beta_pl0.isc_so1(:,3),beta_ot0.isc_ns1(:,3),beta_pl0.isc_ns1(:,3)];
ax1=figureX1([23 20]);
x_pos=figureTmp1_OT(xx1,3);hold on;plot([xlim],[0 0],'k--','LineWidth',2.5);
setFixedPlotArea(ax1, [15 15]);
box off;ylabel('coefficient (∆V)');
set(gca,'XTick',x_pos,'XTickLabel',{'OT','PL','OT','PL'}); 
setAxesForPPT(ax1, 35, 2.5);gcaf1(gca,35);cutBackColor(gca,gcf);
print(gcf,[f_file,'choice_evidence.svg'], '-dsvg','-vector','-r600');
% choice ΣV
xx1=[beta_ot0.isc_so1(:,2),beta_pl0.isc_so1(:,2),beta_ot0.isc_ns1(:,2),beta_pl0.isc_ns1(:,2)];
ax1=figureX1([23 20]);
x_pos=figureTmp1_OT(xx1,3);
hold on;plot([xlim],[0 0],'k--','LineWidth',3);markSignificant(xx1,x_pos);
setFixedPlotArea(ax1, [15 15]);
box off;ylabel('coefficient (∑V)','Interpreter','latex');set(gca,'XTick',x_pos,'XTickLabel',{'OT','PL','OT','PL'}); 
setAxesForPPT(ax1, 35, 2.5);gcaf1(gca,35);cutBackColor(gca,gcf);
print(gcf,[f_file,'choice_bias.svg'], '-dsvg','-vector','-r600');

%% confidence GLM
% confidence ΔV
xx1=[beta_ot0.cf_so20(:,3),beta_pl0.cf_so20(:,3),beta_ot0.cf_ns20(:,3),beta_pl0.cf_ns20(:,3)];
ax1=figureX1([23 20]);
x_pos=figureTmp1_OT(xx1,2);hold on;
setFixedPlotArea(ax1, [15 15]);
box off;ylabel('coefficient (∆V)');set(gca,'XTick',x_pos,'XTickLabel',{'OT','PL','OT','PL'}); gcaf1(gca,35);
setAxesForPPT(ax1, 35, 2.5);cutBackColor(gca,gcf);
print(gcf,[f_file,'confidence_evidence.svg'], '-dsvg','-vector','-r600');
% confidence ΣV
xx1=[beta_ot0.cf_so20(:,2),beta_pl0.cf_so20(:,2),beta_ot0.cf_ns20(:,2),beta_pl0.cf_ns20(:,2)];
ax1=figureX1([23 20]);
x_pos=figureTmp1_OT(xx1,2);
hold on;plot([xlim],[0 0],'k--','LineWidth',3);markSignificant(xx1,x_pos);
setFixedPlotArea(ax1, [15 15]);
box off;ylabel('coefficient (∑V)');set(gca,'XTick',x_pos,'XTickLabel',{'OT','PL','OT','PL'});gcaf1(gca,35);
setAxesForPPT(ax1, 35, 2.5);cutBackColor(gca,gcf);
print(gcf,[f_file,'confidence_bias.svg'], '-dsvg','-vector','-r600');








%%
%%%%%%%%%%%%%%%%%%%%%%%%% function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function markSignificant(xx1,x_pos)
if nargin<2
x_pos=1;
end
ymax_pos=max(xx1,[],'all');
ymin_pos=min(xx1,[],'all');
[~,p_all]=ttest(xx1);
for i=1:size(xx1,2)
    p=p_all(i);
if p < 0.001
    sig_text1 = '***';
elseif p < 0.01
    sig_text1 = '**';
elseif p < 0.05
    sig_text1 = '*';
else
    sig_text1 = 'n.s.';
end
if mean(xx1(:,i))>0
    y_pos=ymax_pos*1.05;
else
    y_pos=ymin_pos*1.05;
end
text(x_pos(i), y_pos , sig_text1, ...
    'HorizontalAlignment', 'center', 'FontSize', 33,'FontName','Arial');%, 'FontWeight', 'bold'
end
end