%%%%%%%%%%%%%%%%%%%%% figure S3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;clc;
% Locate bundled inputs relative to this saved script, independent of pwd.
figureDir = fileparts(mfilename('fullpath'));
codeRoot = fileparts(figureDir);
dataRoot = fullfile(codeRoot,'figure_data');
addpath(fullfile(codeRoot,'utilities'));
f_file = [fullfile(codeRoot,'figure_outputs','figureS3') filesep];
if ~isfolder(f_file)
    mkdir(f_file);
end
ot_color = [226 106 83;242 199 178;88 118 227;171 198 251]./255;
load(fullfile(dataRoot,'behavioral_results','all_GLM_beta_withGroup_fixExpID.mat'),'beta_ot0','beta_pl0');
load(fullfile(dataRoot,'model_results','all_GLM_beta_pred.mat'),'pred_beta');
load(fullfile(dataRoot,'model_results','model_compare_data.mat'),'mc_bic','mc_aic','mc_nll','bmc');

load(fullfile(dataRoot,'model_results','model_parameterRecovery_85.mat'),'ots_minb2','pls_minb2','otn_minb2','pln_minb2');
ots_xx2=cat(1,ots_minb2.optParams);
pls_xx2=cat(1,pls_minb2.optParams);
otn_xx2=cat(1,otn_minb2.optParams);
pln_xx2=cat(1,pln_minb2.optParams);

load(fullfile(dataRoot,'model_results','model_fit_parameter_85.mat'),'ots_minb1','pls_minb1','otn_minb1','pln_minb1');
ots_xx1=cat(1,ots_minb1.optParams);
pls_xx1=cat(1,pls_minb1.optParams);
otn_xx1=cat(1,otn_minb1.optParams);
pln_xx1=cat(1,pln_minb1.optParams);
%% figure S3A
plotModelCompareFigure(mc_aic.os,mc_bic.os,bmc.os,'social-OT',f_file);
plotModelCompareFigure(mc_aic.ps,mc_bic.ps,bmc.ps,'social-PL',f_file);
plotModelCompareFigure(mc_aic.on,mc_bic.on,bmc.on,'nonsocial-OT',f_file);
plotModelCompareFigure(mc_aic.pn,mc_bic.pn,bmc.pn,'nonsocial-PL',f_file);




%% figure S3B
int_color = [132 121 181]./255;int_color1=[132 121 181]./255;
paramNames = {'v_0', 'b_s/b_0', 'b_p/b_0', 'k_s','k_p', 'σ', 'log(λ)', 'tau', 't_0'};
paramNames2 = {'v_0', 'b_1', 'b_2', 'k1','k2', 'σ', 'log(λ)', 'τ', 't_0'};
pn=[7 1:6 8 9];
for i=1:9
    ax1=figureX1([23 21]);
    hold on;
    xx1=[ots_xx1(:,pn(i));pls_xx1(:,pn(i));otn_xx1(:,pn(i));pln_xx1(:,pn(i))];
    xx2=[ots_xx2(:,pn(i));pls_xx2(:,pn(i));otn_xx2(:,pn(i));pln_xx2(:,pn(i))];
    if i==1
        xx1=log(xx1);xx2=log(xx2);
    end
    hold on;
    cond_num=length(xx1)./4;
    for j=1:4
        scatter_tmp(xx1(cond_num*(j-1)+1:cond_num*j),xx2(cond_num*(j-1)+1:cond_num*j),ot_color(j,:));
    end
    x_m=xlim;
    plot([x_m(1),x_m(2)],[x_m(1),x_m(2)],'k--','LineWidth',2);
    pl=errorCIbound(xx1,xx2,int_color1,linspace(x_m(1),x_m(2),100)');pl.LineWidth=3;%,linspace(-0.3,0.2,100)'
    plot([x_m(1),x_m(2)],[x_m(1),x_m(2)],'k--','LineWidth',2);
    textSignificant_all(xx1,xx2);xlim tight;ylim tight;
    
    set(gca, 'LineWidth', 3, 'TickDir', 'in');
    setFixedPlotArea(ax1, [15 15]);
    hold off;box off;ylabel(['recovery (',paramNames{pn(i)},')'],'Interpreter','tex');
    xlabel(['groud truth (',paramNames{pn(i)},')'],'Interpreter','tex');gcaf1(gca,35);
    
    setAxesForPPT(ax1, 35, 2.5);cutBackColor(gca,gcf);
    ax1.YTick = ax1.XTick;
    ax1.YTickLabel = ax1.XTickLabel;
    print(gcf,[f_file,'85_ParameterRecovery-corr(',paramNames2{pn(i)},').svg'], '-dsvg','-vector','-r600');
end


%% figure S3C
% parameter recovery: corr experiment beta, prediction beta
xx2=[pred_beta.os_isc1(:,2);pred_beta.ps_isc1(:,2)];
xx1=[beta_ot0.isc_so1(:,2);beta_pl0.isc_so1(:,2)];
figureScatterCorrExpWithPred(xx1,xx2,'ΣV',ot_color(1:2,:));
print(gcf,[f_file,'85_corr(exp&pred,social choice ΣV).svg'], '-dsvg','-vector','-r600');
%
xx2=[pred_beta.on_isc1(:,2);pred_beta.pn_isc1(:,2)];
xx1=[beta_ot0.isc_ns1(:,2);beta_pl0.isc_ns1(:,2)];
figureScatterCorrExpWithPred(xx1,xx2,'ΣV',ot_color(3:4,:));
print(gcf,[f_file,'85_corr(exp&pred,nonsocial choice ΣV).svg'], '-dsvg','-vector','-r600');
%
xx2=[pred_beta.os_isc1(:,3);pred_beta.ps_isc1(:,3)];
xx1=[beta_ot0.isc_so1(:,3);beta_pl0.isc_so1(:,3)];
figureScatterCorrExpWithPred(xx1,xx2,'ΔV',ot_color(1:2,:));
print(gcf,[f_file,'85_corr(exp&pred,social choice ΔV).svg'], '-dsvg','-vector','-r600');
%
xx2=[pred_beta.on_isc1(:,3);pred_beta.pn_isc1(:,3)];
xx1=[beta_ot0.isc_ns1(:,3);beta_pl0.isc_ns1(:,3)];
figureScatterCorrExpWithPred(xx1,xx2,'ΔV',ot_color(3:4,:));
print(gcf,[f_file,'85_corr(exp&pred,nonsocial choice ΔV).svg'], '-dsvg','-vector','-r600');

%
xx2=[pred_beta.os_conf1(:,2);pred_beta.ps_conf1(:,2)];
xx1=[beta_ot0.cf_so20(:,2);beta_pl0.cf_so20(:,2)];
figureScatterCorrExpWithPred(xx1,xx2,'ΣV',ot_color(1:2,:));
print(gcf,[f_file,'85_corr(exp&pred,social confidence ΣV).svg'], '-dsvg','-vector','-r600');
%
xx2=[pred_beta.on_conf1(:,2);pred_beta.pn_conf1(:,2)];
xx1=[beta_ot0.cf_ns20(:,2);beta_pl0.cf_ns20(:,2)];
figureScatterCorrExpWithPred(xx1,xx2,'ΣV',ot_color(3:4,:));
print(gcf,[f_file,'85_corr(exp&pred,nonsocial confidence ΣV).svg'], '-dsvg','-vector','-r600');
%
xx2=[pred_beta.os_conf1(:,3);pred_beta.ps_conf1(:,3)];
xx1=[beta_ot0.cf_so20(:,3);beta_pl0.cf_so20(:,3)];
figureScatterCorrExpWithPred(xx1,xx2,'ΔV',ot_color(1:2,:));
print(gcf,[f_file,'85_corr(exp&pred,social confidence ΔV).svg'], '-dsvg','-vector','-r600');
%
xx2=[pred_beta.on_conf1(:,3);pred_beta.pn_conf1(:,3)];
xx1=[beta_ot0.cf_ns20(:,3);beta_pl0.cf_ns20(:,3)];
figureScatterCorrExpWithPred(xx1,xx2,'ΔV',ot_color(3:4,:));
print(gcf,[f_file,'85_corr(exp&pred,nonsocial confidence ΔV).svg'], '-dsvg','-vector','-r600');










%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function scatter_tmp(xx1,xx2,ot_color)
hold on
scatter(xx1(1:37),xx2(1:37),'MarkerEdgeColor',[1 1 1],'MarkerFaceColor',ot_color,'SizeData',100);
scatter(xx1(38:56),xx2(38:56),'MarkerEdgeColor',[1 1 1],'MarkerFaceColor',ot_color,'SizeData',105,'Marker','square');
scatter(xx1(57:85),xx2(57:85),'MarkerEdgeColor',[1 1 1],'MarkerFaceColor',ot_color,'SizeData',103,'Marker','^');
end

function textSignificant_all(xx1,xx2)
% [r1(1,1),p1(1,1)]=corr(xx1,xx2);
result = corr_bootstrap(xx1, xx2,'Type', 'Pearson','NBoot', 100000,'Alpha', 0.05,'Seed', 42);
r1 = result.r;p1=result.p_boot;
x_p=max(xlim);y_p=min(ylim);
x_r=range(xlim);y_r=range(ylim);
for i=1:size(r1,1)
    p=p1(i);
if p < 0.001
    sig_text1{i} = '***';
elseif p < 0.01
    sig_text1{i} = '**';
elseif p < 0.05
    sig_text1{i} = '*';
else
    sig_text1{i} = [];
end
end
sig_t3=['{\it r} =',num2str(round(r1,2)),sig_text1{1}];
% sig_tall=sprintf([sig_t3]);
sig_tall=sig_t3;
text(x_p-0.4*x_r, y_p+0.1*y_r , sig_tall, ...
    'HorizontalAlignment', 'left', 'FontSize', 35,'FontName','Arial');
end


function figureScatterCorrExpWithPred(xx1,xx2,label_text,ot_color)
% ot_color = [validatecolor('#367DB0'); validatecolor('#9DC7DD'); ...
%     validatecolor('#3D9F3C'); validatecolor('#9ED17B')];
cond_num=length(xx1)./2;
ax1=figureX1([22 20]);
hold on
for i=1:2
    scatter_tmp(xx1(cond_num*i-cond_num+1:cond_num*i),xx2(cond_num*i-cond_num+1:cond_num*i),ot_color(i,:));
end
[lim_min,lim_max]=getFitXYlim(xx1,xx2);
plot([lim_min, lim_max], [lim_min, lim_max], 'k--','LineWidth',2);
axis square; 
pl=errorCIbound(xx1, xx2, ot_color(1,:), linspace(lim_min, lim_max, 100)');pl.LineWidth=3;
textSignificant_all(xx1,xx2);
setFixedPlotArea(ax1, [15 15]);
hold off;xlabel(['experiment (',label_text,')']);set(gca, 'LineWidth', 3, 'TickDir', 'in'); 
ylabel(['prediction (',label_text,')']);gcaf1(gca,35);setAxesForPPT(ax1, 35, 2.5);
end

function [lim_min,lim_max]=getFitXYlim(xx1,xx2)
% 2. Get global minimum and maximum across all data
all_data = [xx1(:); xx2(:)];
data_min = min(all_data);
data_max = max(all_data);

% ==================== Core optimization: Compact limits and adaptive ticks ====================
% Compute a suitable step from the actual data span(Target approximately 4 ticks, Dividing into 3 intervals)
data_range = data_max - data_min;
rough_step = data_range / 3; 

% Find the nearest"rounded"step (1, 2, 5 sequence), Avoid forcing >= logic, Prevent excessive tick intervals
mag = 10^floor(log10(rough_step));
nice_steps = [1, 2,3,4, 5, 10] * mag;
[~, idx] = min(abs(nice_steps - rough_step));
tick_step = nice_steps(idx);

% [Upper-limit handling]: Continue rounding upward, Ensure a tick at the upper-right maximum label
lim_max = ceil(data_max / tick_step) * tick_step;

% [Lower-limit handling]: **Do not round downward**!Fit closely to the data minimum, Leave only 5% visual margin.
% This avoids forcing the minimum tick and compressing data"into the middle".
lim_min = data_min - (data_range * 0.05);

% [Generate ticks]: Start at the first value greater than or equal to lim_min that lies on the regular grid, Up to lim_max end
first_tick = ceil(lim_min / tick_step) * tick_step;
custom_ticks = round(first_tick : tick_step : lim_max, 6);

% Remove floating-point error
lim_min = round(lim_min, 6);
lim_max = round(lim_max, 6);

% [Adaptive 0 tick]: 
% If 0 Already within the compact axis range, but not sampled by the tick grid, add it 0.
% If 0 Too far from the data range(For example, all data are within 0.6 Up to 1.5), do not add it 0, Avoid concentrating scatter points!
if lim_min <= 0 && lim_max >= 0 && ~any(custom_ticks == 0)
    custom_ticks = sort(unique([custom_ticks, 0]));
end

% 3. Apply the computed coordinate limits
xlim([lim_min, lim_max]);
ylim([lim_min, lim_max]);

% 4. Replace the current axes ticks explicitly
ax = gca;
ax.XTick = custom_ticks;
ax.YTick = custom_ticks;

end

%% model compare
function ots_bic4=plotModelCompareFigure(ots_4m,ots_4m2,bmc_tmp,save_text1,f_file)
%%
m4_color1=[204 164 225; 221 196 237; 236 223 245]./255;
m4_color=[136 73 167; 204 164 225; 221 196 237; 236 223 245]./255;
ots_bic4=ots_4m;

ax1=figureX1([23 20]);
x_tmp=0.6:0.6:1.8;
hold on
for i=2:4  
    xx=(ots_bic4(:,i)-ots_bic4(:,1));
    barh(x_tmp(i-1),mean(xx),'FaceColor','none','BarWidth',0.3,'LineWidth',1.5);

    scatter(xx(1:37),(rand(37,1)-0.5).*0.2+x_tmp(i-1),'MarkerEdgeColor','none','MarkerFaceColor',m4_color1(3,:),'SizeData',100,'Marker','o');
    scatter(xx(38:56),(rand(19,1)-0.5).*0.2+x_tmp(i-1),'MarkerEdgeColor','none','MarkerFaceColor',m4_color1(2,:),'SizeData',105,'marker','square');   
    scatter(xx(57:85),(rand(29,1)-0.5).*0.2+x_tmp(i-1),'MarkerEdgeColor','none','MarkerFaceColor',m4_color1(1,:),'SizeData',103,'marker','^');   

    errorbar(mean(xx),x_tmp(i-1),std(xx)./sqrt(length(xx)),'horizontal', 'k-', 'LineWidth', 3.5,'CapSize',12);
end
plot([0 0],[0.2 2.2],'k','LineWidth',2);
xlabel('AIC difference');set(gca,'YTick',x_tmp,'YTickLabel',{'basic+Urg','basic+MI','basic'});ylim([0.2 2.2]);xlim([-50 150]);
setFixedPlotArea(ax1, [15 15]);gcaf1(gca,35);set(gca, 'LineWidth', 3, 'TickDir', 'in'); setAxesForPPT(ax1, 35, 2.5);
cutBackColor(gca,gcf);print(gcf,[f_file,'85_bar_difference',save_text1,'.svg'], '-dsvg','-vector','-r600');


options=[];
options.DisplayWin=false;
% [~,o_os]=VBA_groupBMC((-ots_4m2./2)',options);
ax1=figureX1([23 20]);
hold on;
bar(bmc_tmp.ep,'FaceColor','none','LineWidth',2);plot([xlim],[0.25 0.25],'k--','LineWidth',2);
hold off;
ylabel('exceedance probabilities');set(gca,'XTick',1:4,'XTickLabel',{'full model','basic+Urg','basic+MI','basic'});
set(gca, 'LineWidth', 3, 'TickDir', 'in'); xlim([0.2 4.7]);ylim([0 1]);
setFixedPlotArea(ax1, [15 15]);box off;gcaf1(gca,35);setAxesForPPT(ax1, 35, 2.5);
cutBackColor(gca,gcf);print(gcf,[f_file,'85_BMC',save_text1,'.svg'], '-dsvg','-vector','-r600');
end
