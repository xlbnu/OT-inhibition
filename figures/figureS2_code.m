%%%%%%%%%%%%%%%%%%%%%%%%%% figure S2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;clc;
% Locate bundled inputs relative to this saved script, independent of pwd.
figureDir = fileparts(mfilename('fullpath'));
codeRoot = fileparts(figureDir);
dataRoot = fullfile(codeRoot,'figure_data');
addpath(fullfile(codeRoot,'utilities'));
f_file = [fullfile(codeRoot,'figure_outputs','figureS2') filesep];
if ~isfolder(f_file)
    mkdir(f_file);
end
load(fullfile(dataRoot,'behavioral_results','all_accauc.mat'),'osub','psub');
load(fullfile(dataRoot,'behavioral_results','all_GLM_beta_withGroup_fixExpID.mat'),'beta_ot0','beta_pl0');
load(fullfile(dataRoot,'behavioral_results','coherence_data.mat'),'coh_data');
%% figure S2A
% perception deicion:coherence
xx1=[coh_data.os_cf,coh_data.ps_cf,coh_data.on_cf,coh_data.pn_cf];
figure4Part(xx1,{'coherence_conflict','coherence'},4,0,f_file);
% perception deicion:auc
xx1=[osub.auc_sc(:,1),psub.auc_sc(:,1),osub.auc_nc(:,1),psub.auc_nc(:,1)];
figure4Part(xx1,{'auc_1st_conflict','auc (1st)'},4,0,f_file);
% perception deicion:rt
xx1=[osub.rts_rdmcf(:,1),psub.rts_rdmcf(:,1),osub.rtn_rdmcf(:,1),psub.rtn_rdmcf(:,1)];
figure4Part(xx1,{'rt_1st_conflict','rt (1st)'},4,0,f_file);
% perception deicion:confidence 
xx1=[osub.conf_cfs(:,1),psub.conf_cfs(:,1),osub.conf_cfn(:,1),psub.conf_cfn(:,1)];
figure4Part(xx1,{'confidence_1st_conflict','confidence(1st)'},4,0,f_file);

%% figure S2B
% social decision:acc 2nd
xx1=[osub.acc_sc(:,2),psub.acc_sc(:,2),osub.acc_nc(:,2),psub.acc_nc(:,2)];
figure4Part(xx1,{'acc_2nd','accuracy (2nd)'},4,0,f_file);
% social decision:auc 2nd
xx1=[osub.auc_sc(:,2),psub.auc_sc(:,2),osub.auc_nc(:,2),psub.auc_nc(:,2)];
figure4Part(xx1,{'auc_2nd','auc (2nd)'},4,0,f_file);
% social decision:congruent trial GLM coefficient(∑V)
xx1=[beta_ot0.cg_so1(:,2),beta_pl0.cg_so1(:,2),beta_ot0.cg_ns1(:,2),beta_pl0.cg_ns1(:,2)];
figure4Part(xx1,{'congruent_∑V','coefficient (∑V)'},4,0,f_file);
% social decision:congruent trial GLM coefficient(∆V)
xx1=[beta_ot0.cg_so1(:,3),beta_pl0.cg_so1(:,3),beta_ot0.cg_ns1(:,3),beta_pl0.cg_ns1(:,3)];
figure4Part(xx1,{'congruent_∆V','coefficient (∆V)'},4,0,f_file);











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
ylim_num=ylim();
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











