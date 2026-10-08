%%%%%%%%%%%%%%%%%%%%%%%%%% figure 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;clc;
% Locate bundled inputs relative to this saved script, independent of pwd.
figureDir = fileparts(mfilename('fullpath'));
codeRoot = fileparts(figureDir);
dataRoot = fullfile(codeRoot,'figure_data');
addpath(fullfile(codeRoot,'utilities'));
f_file = [fullfile(codeRoot,'figure_outputs','figure4') filesep];
if ~isfolder(f_file)
    mkdir(f_file);
end
load(fullfile(dataRoot,'MRS_results','Neurotransmitter_data.mat'),'gaba','glx');% neurotransmitter data
load(fullfile(dataRoot,'behavioral_results','MRS_GLM_beta.mat'),'MRS_ot0','MRS_pl0');

load(fullfile(dataRoot,'model_results','model_fit_parameter_MRS29.mat'),'ots_minb1','pls_minb1','otn_minb1','pln_minb1');
%
ots_xx1=cat(1,ots_minb1.optParams);
pls_xx1=cat(1,pls_minb1.optParams);
otn_xx1=cat(1,otn_minb1.optParams);
pln_xx1=cat(1,pln_minb1.optParams);

%% figure 4E
% confidence bias
xx1=[MRS_ot0.cf_so20(:,2),MRS_pl0.cf_so20(:,2),MRS_ot0.cf_ns20(:,2),MRS_pl0.cf_ns20(:,2)];
ax1=figureX1([23 20]);
x_pos=figureTmp1_OT(xx1,3);
hold on;plot([xlim],[0 0],'k--','LineWidth',2);markSignificant(xx1,x_pos);
setFixedPlotArea(ax1, [15 15]);
box off;ylabel('coefficient (∑V)');set(gca,'XTick',x_pos,'XTickLabel',{'OT','PL','OT','PL'});gcaf1(gca,35);
setAxesForPPT(ax1, 35, 2.5);cutBackColor(gca,gcf);
print(gcf,[f_file,'MRS_confidence_bias.svg'], '-dsvg','-vector','-r600');
%% figure 4F
% parameter lambda
p_ylim=[0;0.2]';
paramNames0 = 'mutual inhibition';
paramNames = 'λ';paramNames2 = 'λ';
pn=[7];
for i=1
    ax1=figureX1([23 20]);
    hold on;
    x_pos=figureTmp1_OT([ots_xx1(:,pn(i)),pls_xx1(:,pn(i)),otn_xx1(:,pn(i)),pln_xx1(:,pn(i))],9); 
    setFixedPlotArea(ax1, [15 15]);
    hold off;box off;set(gca,'XTick',x_pos,'XTickLabel',{'OT','PL','OT','PL'}); gcaf1(gca,35);
    ylabel([paramNames0,' (',paramNames,')'],'FontSize',40,'FontName','Arial');
    setAxesForPPT(ax1, 35, 2.5);
    cutBackColor(gca,gcf);print(gcf,[f_file,'modelParameter(',paramNames2,').svg'], '-dsvg','-vector','-r600');
end
%% figure 4G
% GABA 
xx1=[mean(gaba.ots_cf,2),mean(gaba.pls_cf,2),mean(gaba.otn_cf,2),mean(gaba.pln_cf,2)];
ax1=figureX1([23 20]);
x_pos=figureTmp1_OT(xx1,1);
setFixedPlotArea(ax1, [15 15]);setAxesForPPT(ax1, 35, 3.75)
box off;ylabel('GABA (i.u.)');set(gca,'XTick',x_pos,'XTickLabel',{'OT','PL','OT','PL'}); gcaf1(gca,35);setAxesForPPT(ax1, 35, 2.5);
setAxesForPPT(ax1, 35, 2.5);cutBackColor(gca,gcf);
print(gcf,[f_file,'conflict-GABA.svg'], '-dsvg','-vector','-r600');

%% figure 4H
% corr(confidence bias(ΣV) , lambda) interaction effect
int_color = [132 121 181]./255;int_color1=[132 121 181]./255;
xx1=[MRS_ot0.cf_so20(:,2),MRS_pl0.cf_so20(:,2),MRS_ot0.cf_ns20(:,2),MRS_pl0.cf_ns20(:,2)];
xx2=[xx1(:,1)-xx1(:,3)]-[xx1(:,2)-xx1(:,4)];
lambda_ef=ots_xx1(:,7)-pls_xx1(:,7)-otn_xx1(:,7)+pln_xx1(:,7);
ax1=figureX1([23 20]);
hold on;
scatter(lambda_ef,xx2,'MarkerEdgeColor',[1 1 1],'MarkerFaceColor',int_color,'SizeData',303,'Marker','^');
pl=errorCIbound(lambda_ef,xx2,int_color1,linspace(-0.3,0.15, 100)');pl.LineWidth=3;
hold off;
xlabel('interaction effect (λ)');ylabel(sprintf('interaction effect\n(confidence ΣV)'));gcaf1(gca,35);
setFixedPlotArea(ax1, [15 15]);ylim([-0.4 0.2]);xlim([-0.3 0.15]);xticks(-0.3:0.15:0.15);set(gca, 'LineWidth', 3, 'TickDir', 'in'); 
setAxesForPPT(ax1, 35, 2.5);cutBackColor(gca,gcf);
get2CorrInfor(lambda_ef,xx2);
print(gcf,[f_file,'r(confidence ∑V, lambda)_interactionEffect.svg'], '-dsvg','-vector','-r600');

%% figure 4I
x_gaba=[mean(gaba.ots_cf,2)-mean(gaba.pls_cf,2)-mean(gaba.otn_cf,2)+mean(gaba.pln_cf,2)];
x_v=[MRS_ot0.cf_so20(:,2)-MRS_pl0.cf_so20(:,2)-MRS_ot0.cf_ns20(:,2)+MRS_pl0.cf_ns20(:,2)];
x_lambda=[ots_xx1(:,7)-pls_xx1(:,7)-otn_xx1(:,7)+pln_xx1(:,7)];
mediation_version2((x_gaba), (x_lambda), (x_v), 100000);

%% figure 4J
% corr(GABA , lambda) interaction effect
int_color = [132 121 181]./255;int_color1=[132 121 181]./255;
xx1=ots_xx1(:,7)-otn_xx1(:,7)-(pls_xx1(:,7)-pln_xx1(:,7));
xx2=mean(gaba.ots_cf,2)-mean(gaba.otn_cf,2)-(mean(gaba.pls_cf,2)-mean(gaba.pln_cf,2));
ax1=figureX1([23 22]);
hold on
scatter(xx2,xx1,'MarkerEdgeColor',[1 1 1],'MarkerFaceColor',int_color,'SizeData',303,'Marker','^');
pl=errorCIbound(xx2,xx1,int_color1);pl.LineWidth=3;
setFixedPlotArea(ax1, [15 15]);xlim([-3 3]);
hold off;ylabel('interaction effect (λ)');set(gca, 'LineWidth', 3, 'TickDir', 'in'); xlim([-3 3]);xticks(-3:3);
xlabel(sprintf('interaction effect (GABA) '));gcaf1(gca,35);
setAxesForPPT(ax1, 35, 2.5);cutBackColor(gca,gcf);
get2CorrInfor(xx1,xx2);
print(gcf,[f_file,'r(GABA, lambda)_interactionEffect.svg'], '-dsvg','-vector','-r600');





%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%% function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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



