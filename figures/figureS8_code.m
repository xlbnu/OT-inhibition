%%%%%%%%%%%%%%%%%%%%%%%%% figure S8 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
clear;clc;
% Locate bundled inputs relative to this saved script, independent of pwd.
figureDir = fileparts(mfilename('fullpath'));
codeRoot = fileparts(figureDir);
dataRoot = fullfile(codeRoot,'figure_data');
addpath(fullfile(codeRoot,'utilities'));
f_file = [fullfile(codeRoot,'figure_outputs','figureS8') filesep];
if ~isfolder(f_file)
    mkdir(f_file);
end
ot_color = [226 106 83;242 199 178;88 118 227;171 198 251]./255;
load(fullfile(dataRoot,'MRS_results','Neurotransmitter_data.mat'),'gaba','glx');% neurotransmitter data
load(fullfile(dataRoot,'behavioral_results','MRS_GLM_beta.mat'),'MRS_ot0','MRS_pl0');

load(fullfile(dataRoot,'model_results','model_fit_parameter_MRS29.mat'),'ots_minb1','pls_minb1','otn_minb1','pln_minb1');
ots_xx1=cat(1,ots_minb1.optParams);
pls_xx1=cat(1,pls_minb1.optParams);
otn_xx1=cat(1,otn_minb1.optParams);
pln_xx1=cat(1,pln_minb1.optParams);

lambda_int=[ots_xx1(:,7)-pls_xx1(:,7)-otn_xx1(:,7)+pln_xx1(:,7)];

%% figure S8 C
% phase1 RDM, conflict gaba
xx1=[mean(gaba.ots_cf1,2),mean(gaba.pls_cf1,2),mean(gaba.otn_cf1,2),mean(gaba.pln_cf1,2)];
ax1=figureX1([22 20]);
x_pos=figureTmp1_OT(xx1,1);
setFixedPlotArea(ax1, [15 15]);ylim([3 9]);set(gca,'ytick',3:3:9,'yticklabel',3:3:9);
box off;ylabel('GABA (i.u.)');set(gca,'XTick',x_pos,'XTickLabel',{'OT','PL','OT','PL'}); gcaf1(gca,35);
setAxesForPPT(ax1, 35, 2.5);cutBackColor(gca,gcf);
print(gcf,[f_file,'conflict-GABA-RDM.svg'], '-dsvg','-vector','-r600');
% phase1 RDM, conflict glx
xx1=[mean(glx.ots_cf1,2),mean(glx.pls_cf1,2),mean(glx.otn_cf1,2),mean(glx.pln_cf1,2)];
ax1=figureX1([22 20]);
x_pos=figureTmp1_OT(xx1,1);
setFixedPlotArea(ax1, [15 15]);ylim([10 24]);set(gca,'ytick',10:2:24,'yticklabel',10:2:24);
box off;ylabel('Glx (i.u.)');set(gca,'XTick',x_pos,'XTickLabel',{'OT','PL','OT','PL'}); 
gcaf1(gca,35);setAxesForPPT(ax1, 35, 2.5);cutBackColor(gca,gcf);
print(gcf,[f_file,'conflict-glx-RDM.svg'], '-dsvg','-vector','-r600');
% phase1 RDM, congruent gaba
xx1=[mean(gaba.ots_cg1,2),mean(gaba.pls_cg1,2),mean(gaba.otn_cg1,2),mean(gaba.pln_cg1,2)];
ax1=figureX1([22 20]);
x_pos=figureTmp1_OT(xx1,1);
setFixedPlotArea(ax1, [15 15]);ylim([2 14]);set(gca,'ytick',2:3:14,'yticklabel',2:3:14);
box off;ylabel('GABA (i.u.)');set(gca,'XTick',x_pos,'XTickLabel',{'OT','PL','OT','PL'}); 
gcaf1(gca,35);setAxesForPPT(ax1, 35, 2.5);cutBackColor(gca,gcf);
print(gcf,[f_file,'congruent-GABA-RDM.svg'], '-dsvg','-vector','-r600');

%% figure S8 D
% phase2, conflict gaba
xx1=[mean(gaba.ots_cf,2),mean(gaba.pls_cf,2),mean(gaba.otn_cf,2),mean(gaba.pln_cf,2)];
ax1=figureX1([22 20]);
x_pos=figureTmp1_OT(xx1,1);
setFixedPlotArea(ax1, [15 15]);ylim([3 9]);set(gca,'ytick',3:3:9,'yticklabel',3:3:9);
box off;ylabel('GABA (i.u.)');set(gca,'XTick',x_pos,'XTickLabel',{'OT','PL','OT','PL'}); gcaf1(gca,35);
setAxesForPPT(ax1, 35, 2.5);cutBackColor(gca,gcf);
print(gcf,[f_file,'conflict-GABA.svg'], '-dsvg','-vector','-r600');
% phase2, conflict glx
xx1=[mean(glx.ots_cf,2),mean(glx.pls_cf,2),mean(glx.otn_cf,2),mean(glx.pln_cf,2)];
ax1=figureX1([22 20]);
x_pos=figureTmp1_OT(xx1,1);
setFixedPlotArea(ax1, [15 15]);ylim([10 24]);set(gca,'ytick',10:2:24,'yticklabel',10:2:24);
box off;ylabel('Glx (i.u.)');set(gca,'XTick',x_pos,'XTickLabel',{'OT','PL','OT','PL'}); 
gcaf1(gca,35);setAxesForPPT(ax1, 35, 2.5);cutBackColor(gca,gcf);
print(gcf,[f_file,'conflict-glx-integration.svg'], '-dsvg','-vector','-r600');
% phase 2, congruent gaba
xx1=[mean(gaba.ots_cg,2),mean(gaba.pls_cg,2),mean(gaba.otn_cg,2),mean(gaba.pln_cg,2)];
ax1=figureX1([22 20]);
x_pos=figureTmp1_OT(xx1,1);
setFixedPlotArea(ax1, [15 15]);ylim([2 14]);set(gca,'ytick',2:3:14,'yticklabel',2:3:14);
box off;ylabel(sprintf('GABA (i.u.)'));set(gca,'XTick',x_pos,'XTickLabel',{'OT','PL','OT','PL'}); 
gcaf1(gca,35);setAxesForPPT(ax1, 35, 2.5);cutBackColor(gca,gcf);
print(gcf,[f_file,'congruent-GABA.svg'], '-dsvg','-vector','-r600');

%% figure S8 E
xx1=[mean(gaba.ots_cf-gaba.ots_cg,2),mean(gaba.otn_cf-gaba.otn_cg,2),mean(gaba.pls_cf-gaba.pls_cg,2),mean(gaba.pln_cf-gaba.pln_cg,2)];
ax1=figureX1([22 20]);
x_pos=figureTmp1_OT(xx1,1);
setFixedPlotArea(ax1, [15 15]);yline(0,'k--','LineWidth',2.5);
box off;ylabel('Glx (i.u.)');set(gca,'XTick',x_pos,'XTickLabel',{'OT','PL','OT','PL'}); 
gcaf1(gca,35);setAxesForPPT(ax1, 35, 2.5);cutBackColor(gca,gcf);
print(gcf,[f_file,'conflict-congruent_phase2-GABA.svg'], '-dsvg','-vector','-r600');

%% figure S8 F
x_gaba=[mean(gaba.ots_cf-gaba.ots_cg,2)-mean(gaba.otn_cf-gaba.otn_cg,2)-mean(gaba.pls_cf-gaba.pls_cg,2)+mean(gaba.pln_cf-gaba.pln_cg,2)];
x_v=[MRS_ot0.cf_so20(:,2)-MRS_pl0.cf_so20(:,2)-MRS_ot0.cf_ns20(:,2)+MRS_pl0.cf_ns20(:,2)];
x_lambda=[ots_xx1(:,7)-pls_xx1(:,7)-otn_xx1(:,7)+pln_xx1(:,7)];
mediation_version2((x_gaba), (x_lambda), (x_v), 100000);

%% figure S8 G
% corr(GABA , lambda) interaction effect
int_color = [132 121 181]./255;int_color1=[132 121 181]./255;
xx1=ots_xx1(:,7)-otn_xx1(:,7)-(pls_xx1(:,7)-pln_xx1(:,7));
xx2=mean(gaba.ots_cf-gaba.ots_cg,2)-mean(gaba.otn_cf-gaba.otn_cg,2)-mean(gaba.pls_cf-gaba.pls_cg,2)+mean(gaba.pln_cf-gaba.pln_cg,2);
ax1=figureX1([23 22]);
hold on
scatter(xx2,xx1,'MarkerEdgeColor',[1 1 1],'MarkerFaceColor',int_color,'SizeData',303,'Marker','^');
pl=errorCIbound(xx2,xx1,int_color1);pl.LineWidth=3;
setFixedPlotArea(ax1, [15 15]);xlim([-3 3]);
hold off;ylabel('interaction effect (λ)');set(gca, 'LineWidth', 3, 'TickDir', 'in'); xlim([-3 3]);xticks(-3:3);
xlabel(sprintf('interaction effect (GABA) '));gcaf1(gca,35);
setAxesForPPT(ax1, 35, 2.5);cutBackColor(gca,gcf);
get2CorrInfor(xx2,xx1);
print(gcf,[f_file,'r(GABA, lambda)_interactionEffect.svg'], '-dsvg','-vector','-r600');




