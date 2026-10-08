%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% figure 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Locate bundled inputs relative to this saved script, independent of pwd.
figureDir = fileparts(mfilename('fullpath'));
codeRoot = fileparts(figureDir);
dataRoot = fullfile(codeRoot,'figure_data');
addpath(fullfile(codeRoot,'utilities'));
f_file = [fullfile(codeRoot,'figure_outputs','figure2') filesep];
if ~isfolder(f_file)
    mkdir(f_file);
end
ot_color = [226 106 83;242 199 178; 88 118 227;171 198 251]./255;
load(fullfile(dataRoot,'behavioral_results','all_GLM_beta_withGroup_fixExpID.mat'),'beta_ot0','beta_pl0');
load(fullfile(dataRoot,'model_results','all_GLM_beta_pred.mat'),'pred_beta');

load(fullfile(dataRoot,'model_results','model_fit_parameter_85.mat'),'ots_minb1','pls_minb1','otn_minb1','pln_minb1');
ots_xx1=cat(1,ots_minb1.optParams);
pls_xx1=cat(1,pls_minb1.optParams);
otn_xx1=cat(1,otn_minb1.optParams);
pln_xx1=cat(1,pln_minb1.optParams);

%% figure 2B
% parameter
p_ylim=[-0 0 0 0 0 0 0 0;1.5 1.5 1.5 15 25 0.2 1500 1500]';
paramNames0 = {'baseline drift', 'baseline self bias', 'baseline partner bias', ' self input gain', ' partner input gain', 'noise strength', 'mutual inhibition', 'urgency time constant', 'urgency onset'};
paramNames = {'v_0', 'b_s/b_0', 'b_p/b_0', 'k_s', 'k_p', 'σ', 'λ', 'τ', 't_0'};
paramNames2 = {'v_0', 'b_s', 'b_p', 'k_s','k_p', 'σ', 'λ', 'τ', 't_0'};
pn=[7 1:6 8 9];
for i=1:9
    ax1=figureX1([23 20]);
    hold on;
    x_pos=figureTmp1_OT([ots_xx1(:,pn(i)),pls_xx1(:,pn(i)),otn_xx1(:,pn(i)),pln_xx1(:,pn(i))],9); 
    setFixedPlotArea(ax1, [15 15]);
    
    hold off;box off;set(gca,'XTick',x_pos,'XTickLabel',{'OT','PL','OT','PL'}); gcaf1(gca,35);
    ylabel([paramNames0{pn(i)},' (',paramNames{pn(i)},')'],'FontSize',40,'FontName','Arial');
    setAxesForPPT(ax1, 35, 2.5);cutBackColor(gca,gcf);
    print(gcf,[f_file,'modelParameter(',paramNames2{(pn(i))},').svg'], '-dsvg','-vector','-r600');
end

%% figure 2C
getThreeFactorsInteraction_stat([pred_beta.os_isc1(:,2),pred_beta.ps_isc1(:,2),pred_beta.on_isc1(:,2),pred_beta.pn_isc1(:,2),...
    pred_beta.os_isc1(:,3),pred_beta.ps_isc1(:,3),pred_beta.on_isc1(:,3),pred_beta.pn_isc1(:,3)]);
getThreeFactorsInteraction_stat([pred_beta.os_conf1(:,2),pred_beta.ps_conf1(:,2),pred_beta.on_conf1(:,2),pred_beta.pn_conf1(:,2),...
    pred_beta.os_conf1(:,3),pred_beta.ps_conf1(:,3),pred_beta.on_conf1(:,3),pred_beta.pn_conf1(:,3)]);

% choice ΣV
xx1=[pred_beta.os_isc1(:,2),pred_beta.ps_isc1(:,2),pred_beta.on_isc1(:,2),pred_beta.pn_isc1(:,2)];
ax1=figureX1([23 20]);
x_pos=figureTmp1_OT(xx1,2);plot([xlim],[0 0],'k--','LineWidth',3);markSignificant(xx1,x_pos);
setFixedPlotArea(ax1, [15 15]);
ylabel('coefficient(∑V)');set(gca,'XTick',x_pos,'XTickLabel',{'OT','PL','OT','PL'});gcaf1(gca,35);
setAxesForPPT(ax1, 35, 2.5);cutBackColor(gca,gcf);
print(gcf,[f_file,'model_predict_choice_sigmaV.svg'], '-dsvg','-vector','-r600');
% choice ΔV
xx1=[pred_beta.os_isc1(:,3),pred_beta.ps_isc1(:,3),pred_beta.on_isc1(:,3),pred_beta.pn_isc1(:,3)];
ax1=figureX1([23 20]);
x_pos=figureTmp1_OT(xx1,2);
setFixedPlotArea(ax1, [15 15]);
ylabel('coefficient (∆V)');set(gca,'XTick',x_pos,'XTickLabel',{'OT','PL','OT','PL'});gcaf1(gca,35);
setAxesForPPT(ax1, 35, 2.5);cutBackColor(gca,gcf);
print(gcf,[f_file,'model_predict_choice_deltaV.svg'], '-dsvg','-vector','-r600');
% confidence ΣV
xx1=[pred_beta.os_conf1(:,2),pred_beta.ps_conf1(:,2),pred_beta.on_conf1(:,2),pred_beta.pn_conf1(:,2)];
ax1=figureX1([23 20]);
x_pos=figureTmp1_OT(xx1,2);plot([xlim],[0 0],'k--','LineWidth',3);markSignificant(xx1,x_pos);
setFixedPlotArea(ax1, [15 15]);
ylabel('coefficient(∑V)');set(gca,'XTick',x_pos,'XTickLabel',{'OT','PL','OT','PL'});gcaf1(gca,35);
setAxesForPPT(ax1, 35, 2.5);cutBackColor(gca,gcf);
print(gcf,[f_file,'model_predict_confidence_sigmaV.svg'], '-dsvg','-vector','-r600');
% confidence ΔV
xx1=[pred_beta.os_conf1(:,3),pred_beta.ps_conf1(:,3),pred_beta.on_conf1(:,3),pred_beta.pn_conf1(:,3)];
ax1=figureX1([23 20]);
x_pos=figureTmp1_OT(xx1,2);hold on;plot([xlim],[0 0],'k--','LineWidth',3);
setFixedPlotArea(ax1, [15 15]);
ylabel('coefficient (∆V)');set(gca,'XTick',x_pos,'XTickLabel',{'OT','PL','OT','PL'});gcaf1(gca,35);
setAxesForPPT(ax1, 35, 2.5);cutBackColor(gca,gcf);
print(gcf,[f_file,'model_predict_confidence_deltaV.svg'], '-dsvg','-vector','-r600');

%% figure 2D
%  social： corr (lambda, confidence ΣV)
xx2=[pred_beta.os_conf1(:,2);pred_beta.ps_conf1(:,2)];
xx1=[ots_xx1(:,7);pls_xx1(:,7)];
ax1=figureX1([23 20]);
hold on
for i=1:2
scatter_tmp(xx1(length(ots_xx1)*i-length(ots_xx1)+1:length(ots_xx1)*i),xx2(length(ots_xx1)*i-length(ots_xx1)+1:length(ots_xx1)*i),ot_color(i,:));
end
for i=1:2
pl=errorCIbound(xx1(length(ots_xx1)*i-length(ots_xx1)+1:length(ots_xx1)*i),xx2(length(ots_xx1)*i-length(ots_xx1)+1:length(ots_xx1)*i),ot_color(i,:));
pl.LineWidth=2;
end
setFixedPlotArea(ax1, [15 15],[1.2 0.6]);
hold off;xlabel({'λ'});set(gca, 'LineWidth', 3, 'TickDir', 'in'); 
ylabel(sprintf('model social confidence\n coefficient (∑V)'));gcaf1(gca,35);
get2CorrInfor(xx1(1:85),xx2(1:85));get2CorrInfor(xx1(86:170),xx2(86:170));
setAxesForPPT(ax1, 35, 2.5);cutBackColor(gca,gcf);
print(gcf,[f_file,'corr(lamdba,confidence ∑V) social.svg'], '-dsvg','-vector','-r600');

%   asocial： corr (lambda, confidence ΣV)
xx2=[pred_beta.on_conf1(:,2);pred_beta.pn_conf1(:,2)];
xx1=[otn_xx1(:,7);pln_xx1(:,7)];
ax1=figureX1([23 20]);
hold on
for i=1:2
scatter_tmp(xx1(length(ots_xx1)*i-length(ots_xx1)+1:length(ots_xx1)*i),xx2(length(ots_xx1)*i-length(ots_xx1)+1:length(ots_xx1)*i),ot_color(i+2,:));
end
for i=1:2
pl=errorCIbound(xx1(length(ots_xx1)*i-length(ots_xx1)+1:length(ots_xx1)*i),xx2(length(ots_xx1)*i-length(ots_xx1)+1:length(ots_xx1)*i),ot_color(i+2,:));
pl.LineWidth=2;
end
setFixedPlotArea(ax1, [15 15],[1.2 0.6]);
hold off;xlabel({'λ'});set(gca, 'LineWidth', 3, 'TickDir', 'in'); 
ylabel(sprintf('model asocial confidence\n coefficient (∑V)'));gcaf1(gca,35);
get2CorrInfor(xx1(1:85),xx2(1:85));get2CorrInfor(xx1(86:170),xx2(86:170));
setAxesForPPT(ax1, 35, 2.5);cutBackColor(gca,gcf);
print(gcf,[f_file,'corr(lamdba,confidence ∑V) nonsocial.svg'], '-dsvg','-vector','-r600');

%% figure 2E
% corr(experiment confidence bias(ΣV) , lambda) interaction effect
int_color = [132 121 181]./255;int_color1=[132 121 181]./255;
xx1=[beta_ot0.cf_so20(:,2),beta_pl0.cf_so20(:,2),beta_ot0.cf_ns20(:,2),beta_pl0.cf_ns20(:,2)];
xx2=[xx1(:,1)-xx1(:,3)]-[xx1(:,2)-xx1(:,4)];
lambda_ef=[(ots_xx1(:,7))-(otn_xx1(:,7))]-[(pls_xx1(:,7))-(pln_xx1(:,7))];
ax1=figureX1([23 20]);
hold on;
scatter_tmp(lambda_ef,xx2,int_color);
pl=errorCIbound(lambda_ef,xx2,int_color1,linspace(-0.3,0.2,100)');pl.LineWidth=3;
hold off;
xlabel('interaction effect (λ)');ylabel(sprintf('interaction effect\n(experiment confidence ΣV)'));gcaf1(gca,35);
set(gca, 'LineWidth', 3, 'TickDir', 'in'); 
setFixedPlotArea(ax1, [15 15],[1.2 0.6]);xlim([-0.3 0.2]);xticks(-0.3:0.1:0.2);set(gca,'XTickLabelRotation',0);
get2CorrInfor(lambda_ef,xx2);
setAxesForPPT(ax1, 35, 2.5);cutBackColor(gca,gcf);
print(gcf,[f_file,'r(confidence ∑V, lambda)_interactionEffect.svg'], '-dsvg','-vector','-r600');

% corr(model prediction confidence bias(ΣV) , lambda) interaction effect
int_color = [132 121 181]./255;int_color1=[132 121 181]./255;
xx1=[pred_beta.os_conf1(:,2),pred_beta.ps_conf1(:,2),pred_beta.on_conf1(:,2),pred_beta.pn_conf1(:,2)];
xx2=[xx1(:,1)-xx1(:,3)]-[xx1(:,2)-xx1(:,4)];
lambda_ef=[(ots_xx1(:,7))-(otn_xx1(:,7))]-[(pls_xx1(:,7))-(pln_xx1(:,7))];
ax1=figureX1([23 20]);
hold on;
scatter_tmp(lambda_ef,xx2,int_color);
pl=errorCIbound(lambda_ef,xx2,int_color1,linspace(-0.3,0.2,100)');pl.LineWidth=3;
hold off;
xlabel('interaction effect (λ)');ylabel(sprintf('interaction effect\n(model confidence ΣV)'));gcaf1(gca,35);
set(gca, 'LineWidth', 3, 'TickDir', 'in'); 
setFixedPlotArea(ax1, [15 15],[1.2 0.6]);xlim([-0.3 0.2]);
get2CorrInfor(lambda_ef,xx2);
setAxesForPPT(ax1, 35, 2.5);cutBackColor(gca,gcf);
print(gcf,[f_file,'r(confidence ∑V prediction, lambda)_interactionEffect.svg'], '-dsvg','-vector','-r600');



















%%
%%%%%%%%%%%%%%%%%%%%%%%%%% function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
function ots_meg1=getSameFieldWithHBVdata(ots_bhv,ots_meg)
ots_field1=fields(ots_bhv);
ots_field2=fields(ots_meg);
ots_field=intersect(ots_field1, ots_field2);
ots_meg1=[];
for i=1:length(ots_meg)
for j=1:length(ots_field)
ots_meg_tmp.(string(ots_field(j)))=ots_meg(i).(string(ots_field(j)));
end
ots_meg1=cat(1,ots_meg1,ots_meg_tmp);
end
end

%
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
%
function scatter_tmp(xx1,xx2,ot_color)
hold on
scatter(xx1(1:37),xx2(1:37),'MarkerEdgeColor',ot_color,'MarkerFaceColor',ot_color,'SizeData',100);
scatter(xx1(38:56),xx2(38:56),'MarkerEdgeColor',ot_color,'MarkerFaceColor',ot_color,'SizeData',105,'Marker','square');
scatter(xx1(57:85),xx2(57:85),'MarkerEdgeColor',ot_color,'MarkerFaceColor',ot_color,'SizeData',103,'Marker','^');
end
