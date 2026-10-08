%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% figure S1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;clc;
% Locate bundled inputs relative to this saved script, independent of pwd.
figureDir = fileparts(mfilename('fullpath'));
codeRoot = fileparts(figureDir);
dataRoot = fullfile(codeRoot,'figure_data');
addpath(fullfile(codeRoot,'utilities'));
f_file = [fullfile(codeRoot,'figure_outputs','figureS1') filesep];
if ~isfolder(f_file)
    mkdir(f_file);
end
load(fullfile(dataRoot,'behavioral_results','all_accauc.mat'),'osub','psub');
load(fullfile(dataRoot,'behavioral_results','singleExp_GLM_beta.mat'),'sbeta_ot0','sbeta_pl0');


%% figure S1A
% synergy effect (2HBT1)
osc=max([osub.acc_sc(:,1),osub.acc_oc(:,1)],[],2);psc=max([psub.acc_sc(:,1),psub.acc_oc(:,1)],[],2);
onc=max([osub.acc_nc(:,1),osub.acc_ic(:,1)],[],2);pnc=max([psub.acc_nc(:,1),psub.acc_ic(:,1)],[],2);
xx1=[osub.acc_sc(:,2)-osc,psub.acc_sc(:,2)-psc,osub.acc_nc(:,2)-onc,psub.acc_nc(:,2)-pnc];
figure4Part(xx1,{'2HBT1','synergy effect'},1,0,f_file);

%% figure S1B
% choice GLM coefficent of each cohort
xx1=[sbeta_ot0.isc_so1(:,3),sbeta_pl0.isc_so1(:,3),sbeta_ot0.isc_ns1(:,3),sbeta_pl0.isc_ns1(:,3)];
figure4Part(xx1,{'choice_∆V','coefficient (∆V)'},1,0,f_file);
xx1=[sbeta_ot0.isc_so1(:,2),sbeta_pl0.isc_so1(:,2),sbeta_ot0.isc_ns1(:,2),sbeta_pl0.isc_ns1(:,2)];
figure4Part(xx1,{'choice_∑V','coefficient (∑V)'},1,0,f_file);

%% figure S1C
% confidence GLM coefficent of each cohort
xx1=[sbeta_ot0.cf_so20(:,3),sbeta_pl0.cf_so20(:,3),sbeta_ot0.cf_ns20(:,3),sbeta_pl0.cf_ns20(:,3)];
figure4Part(xx1,{'confidence_∆V','coefficient (∆V)'},1,0,f_file);
xx1=[sbeta_ot0.cf_so20(:,2),sbeta_pl0.cf_so20(:,2),sbeta_ot0.cf_ns20(:,2),sbeta_pl0.cf_ns20(:,2)];
figure4Part(xx1,{'confidence_∑V','coefficient (∑V)'},1,0,f_file);












%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function figure4Part(xx1,save_text,x_group,n_chance,f_file)
if nargin<3
x_group=1;
end
if nargin<4
n_chance=0;
end
save_text1=save_text{1};
ylabel_text=save_text{2};
exp_idx=[ones(37,1);ones(19,1).*2;ones(29,1).*3];

% BHV sub37
xx1_tmp=xx1(exp_idx==1,:);
ax1=figureX1([22 20]);
x_pos=figureTmp1_OT(xx1_tmp,x_group);
ylim_num=ylim;
if ylim_num(1)<0 && ylim_num(2)>0 
plot(xlim,[n_chance n_chance],'k--','LineWidth',2);
end
setFixedPlotArea(ax1, [15 15]);
box off;ylabel(ylabel_text,'Interpreter','tex');set(gca,'XTick',x_pos,'XTickLabel',{'OT','PL','OT','PL'});gcaf1(gca,35);
setAxesForPPT(ax1, 35, 3.75);cutBackColor(gca,gcf);
print(gcf,[f_file,'37_',save_text1,'.svg'], '-dsvg','-vector','-r600');
% MEG sub19
xx1_tmp=xx1(exp_idx==2,:);
ax1=figureX1([22 20]);
x_pos=figureTmp1_OT(xx1_tmp,x_group);
if ylim_num(1)<0 && ylim_num(2)>0 
plot(xlim,[n_chance n_chance],'k--','LineWidth',2);
end
setFixedPlotArea(ax1, [15 15]);
box off;ylabel(ylabel_text,'Interpreter','tex');set(gca,'XTick',x_pos,'XTickLabel',{'OT','PL','OT','PL'});gcaf1(gca,35);
setAxesForPPT(ax1, 35, 3.75);cutBackColor(gca,gcf);
print(gcf,[f_file,'19_',save_text1,'.svg'], '-dsvg','-vector','-r600');
% fMRS sub29
xx1_tmp=xx1(exp_idx==3,:);
ax1=figureX1([22 20]);
x_pos=figureTmp1_OT(xx1_tmp,x_group);
if ylim_num(1)<0 && ylim_num(2)>0 
plot(xlim,[n_chance n_chance],'k--','LineWidth',2);
end
setFixedPlotArea(ax1, [15 15]);
box off;ylabel(ylabel_text,'Interpreter','tex');set(gca,'XTick',x_pos,'XTickLabel',{'OT','PL','OT','PL'});gcaf1(gca,35);
setAxesForPPT(ax1, 35, 3.75);cutBackColor(gca,gcf);
print(gcf,[f_file,'29_',save_text1,'.svg'], '-dsvg','-vector','-r600');
end
