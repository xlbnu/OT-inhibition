
%%%%%%%%%%%%%%%%%%%%%%%%%%%% figure 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;clc;
% Locate bundled inputs relative to this saved script, independent of pwd.
figureDir = fileparts(mfilename('fullpath'));
codeRoot = fileparts(figureDir);
dataRoot = fullfile(codeRoot,'figure_data');
addpath(fullfile(codeRoot,'utilities'));
f_file = [fullfile(codeRoot,'figure_outputs','figure3') filesep];
if ~isfolder(f_file)
    mkdir(f_file);
end
assert(exist('ft_freqgrandaverage','file') == 2, ...
    'Configure FieldTrip first: run setup_figure_paths(fieldtripRoot) from the package root.');
ot_color = [226 106 83;242 199 178;88 118 227;171 198 251]./255;
% bhv 
load(fullfile(dataRoot,'behavioral_results','MEG_GLM_beta.mat'),'MEG_ot0','MEG_pl0');
load(fullfile(dataRoot,'behavioral_results','MEG_acc_auc.mat'),'osub','psub');
% model
load(fullfile(dataRoot,'model_results','model_fit_parameter_MEG19.mat'));
ots_xx1=cat(1,ots_minb1.optParams);
pls_xx1=cat(1,pls_minb1.optParams);
otn_xx1=cat(1,otn_minb1.optParams);
pln_xx1=cat(1,pln_minb1.optParams);




%% figure3B
% synergy effect (2HBT1)
osc=max([osub.acc_sc(:,1),osub.acc_oc(:,1)],[],2);psc=max([psub.acc_sc(:,1),psub.acc_oc(:,1)],[],2);
onc=max([osub.acc_nc(:,1),osub.acc_ic(:,1)],[],2);pnc=max([psub.acc_nc(:,1),psub.acc_ic(:,1)],[],2);
xx1=[osub.acc_sc(:,2)-osc,psub.acc_sc(:,2)-psc,osub.acc_nc(:,2)-onc,psub.acc_nc(:,2)-pnc];
ax1=figureX1([23 20]);
x_pos=figureTmp1_OT(xx1,3);hold on;plot([xlim],[0 0],'k--','LineWidth',2);
box off;ylabel('synergy effect');set(gca,'XTick',x_pos,'XTickLabel',{'OT','PL','OT','PL'}); gcaf1(gca,35);
setFixedPlotArea(ax1, [15 15]);setAxesForPPT(ax1, 35, 2.5);cutBackColor(gca,gcf);
print(gcf,[f_file,'acc_2HBT1.svg'], '-dsvg','-vector','-r600');

%% figure3C
% confidence GLM beta(ΣV)(confirmation bias)
xx1=[MEG_ot0.cf_so20(:,2),MEG_pl0.cf_so20(:,2),MEG_ot0.cf_ns20(:,2),MEG_pl0.cf_ns20(:,2)];
ax1=figureX1([24 20]);
x_pos=figureTmp1_OT(xx1,1);
hold on;plot([xlim],[0 0],'k--');markSignificant(xx1,x_pos);
setFixedPlotArea(ax1, [15 15]);
box off;ylabel('coefficient (∑V)');set(gca,'XTick',x_pos,'XTickLabel',{'OT','PL','OT','PL'});gcaf1(gca,35);
setAxesForPPT(ax1, 35, 2.5);cutBackColor(gca,gcf);
print(gcf,[f_file,'confidence_bias.svg'], '-dsvg','-vector','-r600');

%% figure3D
% model parameter
paramNames0 = 'mutual inhibition';
paramNames = 'λ';
pn=[7];
for i=1
    ax1=figureX1([24 20]);
    hold on;
    x_pos=figureTmp1_OT([ots_xx1(:,pn(i)),pls_xx1(:,pn(i)),otn_xx1(:,pn(i)),pln_xx1(:,pn(i))],9); 
    setFixedPlotArea(ax1, [15 15]);    
    hold off;box off;set(gca,'XTick',x_pos,'XTickLabel',{'OT','PL','OT','PL'}); gcaf1(gca,35);
    ylabel([paramNames0,' (',paramNames,')'],'FontSize',40,'FontName','Arial');
    setAxesForPPT(ax1, 35, 2.5);cutBackColor(gca,gcf);
    print(gcf,[f_file,'modelParameter(',paramNames,').svg'], '-dsvg','-vector','-r600');
end

%% figure3E
% topoplot ΣV
load(fullfile(dataRoot,'MEG_results','pressPutton_task_and_ITI_glm_beta_OT_PL.mat'));
tidx_base = efr_glm.baseline_time >= 0 & efr_glm.baseline_time <= 0.2;
base_iti=squeeze(mean(efr_glm.beta_bias_base(:,:,:,tidx_base,2),[1 4],'omitnan'));
tidx = efr_glm.sensor_time >= -0.5 & efr_glm.sensor_time <= -0.3;
sub_beta = squeeze(mean(efr_glm.beta_bias(:,:,:,tidx,2), [1 4], 'omitnan'))-base_iti;
% Wrap in FieldTrip timelock format
topo = [];
topo.label = regexprep(efr_glm.sensor_label(:), '\+MEG', '+');
topo.time = 0;
topo.avg = mean(sub_beta, 1, 'omitnan')';
topo.avg_base = mean(base_iti, 1, 'omitnan')';
topo.dimord = 'chan_time';
% 
cfg = [];
cfg.xlim         =[-0.5 -0.3];
cfg.zlim         = [-0.02 0.02];
cfg.marker       = 'off';
cfg.layout       = 'neuromag306cmb.lay';
cfg.comment      = 'no';
cfg.colorbar     = 'no';
cfg.parameter = 'avg';
cmap = get_adaptive_cmap(cfg.zlim);
cfg.figure=figureX1([23 20]);
ft_topoplotER(cfg, topo);
colormap(cmap);
setFixedPlotArea(cfg.figure, [15 15]);
print(gcf,[f_file,'sensorGLM_coefficent(ΣV)_topgraphic'], '-dsvg','-r600');

%% figure3F
% brainNet to generate coefficient(ΣV) source ∩ oxytocin reciptor gene
% expression map
cfg = struct();
cfg.brainnet_dir = fullfile(codeRoot,'utilities','BrainNet-Viewer');
cfg.surface_file = fullfile(cfg.brainnet_dir,'BrainMesh_ICBM152Right.nv');
cfg.cfg_file = fullfile(dataRoot,'MEG_results','BrainNet_cfg_medial_right.mat');

cfg.mapping_dir = fullfile(dataRoot,'MEG_results');
cfg.output_dir = f_file;

% Render the map assigned to this panel.
cfg.input_pattern = 'decision_sigmaV_source_OXTRconjunction_.nii';

% Preserve saved display thresholds and color limits.
cfg.fix_global_color_scale = false;
cfg.save_tif_only = true;
result = batch_brainnet_nii_to_tif(cfg);

%% figure3G-I
load(fullfile(dataRoot,'model_results','model_DDM_trajectory_GLM.mat'),'beta_acc','beta_acc1','bs_acc','bc_acc','bs_acc1','bc_acc1');
%% figure3G 
% trajectory
trajectory_rang = 10:10:1000;
staySelf  = smoothdata(squeeze(mean(bs_acc(:,:,trajectory_rang,1),1,'omitnan')),2,'gaussian',50);
stayOther = smoothdata(squeeze(mean(bs_acc(:,:,trajectory_rang,2),1,'omitnan')),2,'gaussian',50);
changeSelf  = smoothdata(squeeze(mean(bc_acc(:,:,trajectory_rang,1),1,'omitnan')),2,'gaussian',50);
changeOther = smoothdata(squeeze(mean(bc_acc(:,:,trajectory_rang,2),1,'omitnan')),2,'gaussian',50);

staySelf1  = smoothdata(squeeze(mean(bs_acc1(:,:,trajectory_rang,1),1,'omitnan')),2,'gaussian',50);
stayOther1 = smoothdata(squeeze(mean(bs_acc1(:,:,trajectory_rang,2),1,'omitnan')),2,'gaussian',50);
changeSelf1  = smoothdata(squeeze(mean(bc_acc1(:,:,trajectory_rang,1),1,'omitnan')),2,'gaussian',50);
changeOther1 = smoothdata(squeeze(mean(bc_acc1(:,:,trajectory_rang,2),1,'omitnan')),2,'gaussian',50);

lightColor = [187 181 214]./255;
darkColor  = [132 121 181]./255;
grayLight  = [0.82 0.82 0.82];
grayDark   = [0.25 0.25 0.25];

ax1 = figureX1([24 21]);
hold(ax1,'on');

% without lambda: stay
hStayNoLambda = plot_gradient_trajectory_sem( ...
    staySelf1,stayOther1,grayLight,grayDark, ...
    semMultiplier=1, ...
    semAlpha=0.18, ...
    semColorStrength=0.45, ...
    lineWidth=8, ...
    startMarker='o', ...
    endMarker='o', ...
    startMarkerSize=280, ...
    endMarkerSize=280, ...
    markerLineWidth=1.5);

% without lambda: change
hChangeNoLambda = plot_gradient_trajectory_sem( ...
    changeSelf1,changeOther1,grayLight,grayDark, ...
    semMultiplier=1,semAlpha=0.18, ...
    lineWidth=8, ...
    startMarker='o',endMarker='o', ...
    startMarkerSize=280,endMarkerSize=280);

% with lambda: stay
hStay = plot_gradient_trajectory_sem( ...
    staySelf,stayOther,lightColor,darkColor, ...
    semMultiplier=1,semAlpha=0.25, ...
    lineWidth=8, ...
    startMarker='o',endMarker='o', ...
    startMarkerSize=280,endMarkerSize=280);

% with lambda: change
hChange = plot_gradient_trajectory_sem( ...
    changeSelf,changeOther,lightColor,darkColor, ...
    semMultiplier=1,semAlpha=0.25, ...
    lineWidth=8, ...
    startMarker='o',endMarker='o', ...
    startMarkerSize=280,endMarkerSize=280);

xlabel(ax1,'magnitude(self)');
ylabel(ax1,'magnitude(partner)');
xlim(ax1,[0 0.08]);
ylim(ax1,[0 0.08]);

hline = refline(ax1,1,0);
hline.LineStyle = '--';
hline.Color = 'k';
hline.LineWidth = 1.5;

xline(ax1,0,'k--');
yline(ax1,0,'k--');
axis(ax1,'square');
box(ax1,'off');

setFixedPlotArea(ax1,[15 15]);
setAxesForPPT(ax1,35,2.5);
gcaf1(ax1,35);

print(gcf,[f_file,'model_trajectory_self-other.svg'], '-dsvg','-vector','-r600');

%% figure3H
% social/asocial : OT(coefficient (ΣV)) - PLcoefficient (ΣV)   
int_color1=[132 121 181]./255;
ax1=figureX1([25 21]);
differ_social=squeeze(beta_acc(1,:,:,2)-beta_acc(2,:,:,2));differ_asocial=squeeze(beta_acc(3,:,:,2)-beta_acc(4,:,:,2));
differ_social1=squeeze(beta_acc1(1,:,:,2)-beta_acc1(2,:,:,2));differ_asocial1=squeeze(beta_acc1(3,:,:,2)-beta_acc1(4,:,:,2));
hold on
errorbound(mean(differ_social1,1,'omitnan'),std(differ_social1,[],1,'omitnan')./sqrt(sum(~isnan(differ_social1))),'color',[145 91 78]./255,'linewidth',5);%
errorbound(mean(differ_asocial1,1,'omitnan'),std(differ_asocial1,[],1,'omitnan')./sqrt(sum(~isnan(differ_asocial1))),'color',[72 86 135]./255,'linewidth',5);%
errorbound(mean(differ_social,1,'omitnan'),std(differ_social,[],1,'omitnan')./sqrt(sum(~isnan(differ_social))),'color',ot_color(1,:),'linewidth',5);
errorbound(mean(differ_asocial,1,'omitnan'),std(differ_asocial,[],1,'omitnan')./sqrt(sum(~isnan(differ_asocial))),'color',ot_color(3,:),'linewidth',5);
plot_column_significance(1:1000,differ_social-differ_asocial,'Color',int_color1,'YOffsetRatio',0.20,'LineWidth',5);
plot_column_significance(1:1000,differ_social1-differ_asocial1,'Color',[0.7 0.7 0.7],'YOffsetRatio',-0.20,'LineWidth',5);
yline(0,'k--','LineWidth',1.5);xlim([200 1000]);
hLine = findobj(gca, 'Type', 'Line');hLine(end-5).Color(4) = 0.3;hLine(end-2).Color(4) = 0.3;
set(gca,'xtick',200:200:1000,'XTickLabel',-1000+200:200:1000,'XTickLabelRotation',0);
hold off;xlabel('time from choice made');ylabel(['coefficient (ΣV)',newline,'(OT — PL)']);ylim([-0.6 0.4]);yticks(-0.6 :0.2:0.4);
setFixedPlotArea(ax1, [15 15]);
setAxesForPPT(ax1, 35, 2.5);gcaf1(gca,35);
print(gcf,[f_file,'MEG-BHV-ΣV-interaction-timeCourse.svg'], '-dsvg','-vector','-r600');

%% figure3I
% end point effect
end_wl0=nan(19,4);
for c=1:4
end_wl0(:,c) = squeeze(beta_acc(c,:,1000,2))';
end
%
xx1=end_wl0;
ax1=figureX1([23 20]);
x_pos=figureTmp1_OT(xx1,1);hold on;plot([xlim],[0 0],'k--','LineWidth',3);
setFixedPlotArea(ax1, [15 15]);
box off;ylabel(['coefficient (ΣV)']);
set(gca,'XTick',x_pos,'XTickLabel',{'OT','PL','OT','PL'}); 
setAxesForPPT(ax1, 35, 2.5);gcaf1(gca,35);
print(gcf,[f_file,'MEG-BHV-ΣV-interaction-endPoint.svg'], '-dsvg','-vector','-r600');

%% figure3J-L
load(fullfile(dataRoot,'MEG_results','trajectory_sourceEFR_GLM_beta_sumq.mat'),'sb_chosen0','sb_chosen','cb_chosen0','cb_chosen',...
    'otsbs_chosen','otsbc_chosen','mb_chosen','otsb_chosen1','plsb_chosen1','otnb_chosen1','plnb_chosen1','vchan_time');
%% figure3J
int_color1 = [187 181 214]./255;
int_color  = [132 121 181]./255;

grayLight = [224 218 208]./255;
grayDark  = [ 78  69  59]./255;

stay_sem_color   = [0.88 0.85 0.94];
change_sem_color = [0.88 0.85 0.94];
mixed_sem_color  = [0.86 0.86 0.86];

% =========================================================
% Single-participant signed beta: subject × time
% ==========================================================
stay_self_sub = squeeze(otsbs_chosen(:,:,2));
stay_other_sub = squeeze(otsbs_chosen(:,:,3));

change_self_sub = squeeze(otsbc_chosen(:,:,2));
change_other_sub = squeeze(otsbc_chosen(:,:,3));

mixed_self_sub = squeeze(mb_chosen(:,:,2));
mixed_other_sub = squeeze(mb_chosen(:,:,3));

% time window
stay_idx = vchan_time>-0.7 & vchan_time<-0.35;
change_idx = vchan_time>-0.4 & vchan_time<-0.10;
mixed_idx = vchan_time>-0.7 & vchan_time<-0.0;

% =========================================================
% Bootstrap trajectory
% ==========================================================
stay_result = bootstrap_trajectory_sem(stay_self_sub,stay_other_sub,stay_idx, ...
    'NBoot',5000, 'FirstSmoothMethod','movmean', 'FirstSmoothWindow',20, ...
    'SecondSmoothMethod','gaussian','SecondSmoothWindow',20,'Seed',2026);

change_result = bootstrap_trajectory_sem(change_self_sub,change_other_sub,change_idx, ...
    'NBoot',5000,'FirstSmoothMethod','movmean','FirstSmoothWindow',20, ...
    'SecondSmoothMethod','movmean','SecondSmoothWindow',20,'Seed',2027);

mixed_result = bootstrap_trajectory_sem(mixed_self_sub,mixed_other_sub,mixed_idx, ...
    'NBoot',5000,'FirstSmoothMethod','movmean','FirstSmoothWindow',20, ...
    'SecondSmoothMethod','movmean','SecondSmoothWindow',20,'Seed',2028);

% =========================================================
% Generate gradient colors correctly
% ==========================================================
fraction_s = linspace(0,1,nnz(stay_idx))';
gradientColors_s = int_color1 + fraction_s.*(int_color-int_color1);

fraction_c = linspace(0,1,nnz(change_idx))';
gradientColors_c = int_color1 + fraction_c.*(int_color-int_color1);

fraction_m = linspace(0,1,nnz(mixed_idx))';
gradientColors_m = grayLight + fraction_m.*(grayDark-grayLight);

% =========================================================
% Plotting
% ==========================================================
ax1 = figureX1([24 21]);
hold(ax1,'on');

% First draw the SEM region, then the mean line
draw_trajectory_sem_union(ax1, mixed_result.upper_x,mixed_result.upper_y, ...
    mixed_result.lower_x,mixed_result.lower_y,mixed_sem_color,0.8,true);

draw_trajectory_sem_union(ax1,stay_result.upper_x,stay_result.upper_y, ...
    stay_result.lower_x,stay_result.lower_y,stay_sem_color,0.35,true);

draw_trajectory_sem_union(ax1,change_result.upper_x,change_result.upper_y, ...
    change_result.lower_x,change_result.lower_y,change_sem_color,0.35,true);

% Gradient-colored mean trajectory
plot_gradient_path(ax1,mixed_result.mean_x,mixed_result.mean_y,gradientColors_m,false,8);

plot_gradient_path(ax1,stay_result.mean_x,stay_result.mean_y,gradientColors_s,false,8);

plot_gradient_path(ax1,change_result.mean_x,change_result.mean_y,gradientColors_c,false,8);

xlabel(ax1,'magnitude(self)');
ylabel(ax1,'magnitude(partner)');

xline(ax1,0,'k--');
yline(ax1,0,'k--');

ylim(ax1,[0 0.06]);
xlim(ax1,[0 0.06]);

% self = other
plot(ax1,[0 0.06],[0 0.06], ...
    'k--','LineWidth',1.5, ...
    'HandleVisibility','off');

yticks(ax1,[0:0.02:0.06]);
xticks(ax1,[0:0.02:0.06]);

axis(ax1,'square');
box(ax1,'off');

setFixedPlotArea(ax1,[15 15]);gcaf1(ax1,35);
setAxesForPPT(ax1,36,2);
hold(ax1,'off');
print(gcf,[f_file,'MEG-source-conjunction-trajectory_sem.svg'], '-dsvg','-vector','-r600');
%% figure3K
% interaction time course
ot_color1 = [226 106 83;242 199 178; 88 118 227;171 198 251]./255;int_color1=[132 121 181]./255;
diff_stay_social = smoothdata(squeeze((otsb_chosen1(:,:,2))),2,'movmean',20)-smoothdata(squeeze((plsb_chosen1(:,:,2))),2,'movmean',20);
diff_stay_asocial = smoothdata(squeeze((otnb_chosen1(:,:,2))),2,'movmean',20)-smoothdata(squeeze((plnb_chosen1(:,:,2))),2,'movmean',20);
ax1=figureX1([25 21]);
hold on
errorbound(vchan_time,mean(diff_stay_social,1,'omitnan'),std(diff_stay_social,'omitnan')./sqrt(size(diff_stay_social,1)),'color',ot_color1(1,:),'linewidth',5);
errorbound(vchan_time,mean(diff_stay_asocial,1,'omitnan'),std(diff_stay_asocial,'omitnan')./sqrt(size(diff_stay_asocial,1)),'color',ot_color1(3,:),'linewidth',5);
plot_column_significance(vchan_time,diff_stay_social-diff_stay_asocial,'Color',int_color1,'YOffsetRatio',0.4,'LineWidth',5,'TimeWindow',[-1 0.5]);
yline(0,'k--','LineWidth',1.5);xlim([-1 0.5]);xline(0,'k--','LineWidth',1.5);%ylim([-0.1 0.1]);
hold off;xlabel('time from choice made');ylabel(['coefficient (ΣV)',newline,'(OT — PL)']);
setFixedPlotArea(ax1, [15 15]);
gcaf1(gca,35);setAxesForPPT(ax1, 36, 2.5);cutBackColor(gca,gcf);
print(gcf,[f_file,'MEG-source-conjunction-interaction-timeCourse.svg'], '-dsvg','-vector','-r600');
%% figure3L
% interaction
time_idx0=vchan_time>-0.4 & vchan_time<-0.3;
xx2=[mean(smoothdata(squeeze(otsb_chosen1(:,time_idx0,2)),2,'movmean',20),2),...
    mean(smoothdata(squeeze(plsb_chosen1(:,time_idx0,2)),2,'movmean',20),2),...
    mean(smoothdata(squeeze(otnb_chosen1(:,time_idx0,2)),2,'movmean',20),2),...
    mean(smoothdata(squeeze(plnb_chosen1(:,time_idx0,2)),2,'movmean',20),2)];
ax1=figureX1([23 20]);
x_pos=figureTmp1_OT(xx2,1);hold on;plot([xlim],[0 0],'k--','LineWidth',2.5);
setFixedPlotArea(ax1, [15 15]);
box off;ylabel('coefficient (ΣV)');
set(gca,'XTick',x_pos,'XTickLabel',{'OT','PL','OT','PL'}); 
gcaf1(gca,35);setAxesForPPT(ax1, 36, 2.5);cutBackColor(gca,gcf);
print(gcf,[f_file,'MEG-source-conjunction-interaction_plot.svg'], '-dsvg','-vector','-r600');


%% figure3 M-P
% mean-field model load data
load(fullfile(dataRoot,'MEG_results','40_level_data.mat'))
load(fullfile(dataRoot,'MEG_results','pair_data.mat'),'order_idx1','order_idx2');
low_mtfr=squeeze(data(22:40,:,:,:));inhibit_low=inhibition_strength(22:40);
high_mtfr=squeeze(data(1:19,:,:,:));inhibit_high=inhibition_strength(1:19);
avg_mtfr=(mean((low_mtfr+high_mtfr)./2,1));
diff_mtfr=(mean(low_mtfr-high_mtfr,1));
f_adj=f;

mtfr_avg=[];
mtfr_avg.time=double(t_wind_start+100)./1000;
mtfr_avg.freq=f_adj;
mtfr_avg.powspctrm=avg_mtfr;
mtfr_avg.dimord='chan_freq_time';
mtfr_avg.label={'cluster1'};

mtfr_diff=[];
mtfr_diff.time=double(t_wind_start+100)./1000;
mtfr_diff.freq=f_adj;
mtfr_diff.powspctrm=diff_mtfr;
mtfr_diff.dimord='chan_freq_time';
mtfr_diff.label={'cluster1'};

%% figure3M 
avg_smooth=getSmoothTFR(mtfr_avg);
cfg_plot_diff = [];
cfg_plot_diff.parameter='powspctrm';
cfg_plot_diff.title        =' ';
cfg_plot_diff.xlim=[-1 0.5];
cfg_plot_diff.ylim=[30 120];
cfg_plot_diff.zlim = [0 0.9];
cfg_plot_diff.colorbar     ='no';
cmap=get_adaptive_cmap(cfg_plot_diff.zlim,256,1);
[ax1,cfg_plot_diff.figure]=figureX1([23 20]);
ft_singleplotTFR(cfg_plot_diff, avg_smooth);
colormap(cmap);
hold on;xline(0,'k--','LineWidth',2);
setFixedPlotArea(ax1, [15 15]);set(gca, 'LineWidth', 2); xticks(-1:0.5:0.5);
xlabel('time from choice made');ylabel(sprintf('frequency'));gcaf1(gca,35);
setAxesForPPT(gca, 35, 2.5);cutBackColor(gca,gcf);
print(gcf,[f_file,'modelGamma_avg_weak_strong.svg'], '-dsvg','-vector','-r600');
%% figure3N
diff_smooth=getSmoothTFR(mtfr_diff);
cfg_plot_diff = [];
cfg_plot_diff.parameter='powspctrm';
cfg_plot_diff.title        =' ';
cfg_plot_diff.xlim=[-1 0.5];
cfg_plot_diff.ylim=[30 120];
cfg_plot_diff.zlim = [-0.4 0.4];
cfg_plot_diff.colorbar     ='no';
cmap=get_adaptive_cmap(cfg_plot_diff.zlim);
[ax1,cfg_plot_diff.figure]=figureX1([23 20]);
ft_singleplotTFR(cfg_plot_diff, diff_smooth);
colormap(cmap);
hold on;plot([0 0],[30 120],'k--','LineWidth',2);
setFixedPlotArea(ax1, [15 15]);set(gca, 'LineWidth', 2); xticks(-1:0.5:0.5);
xlabel('time from choice made');ylabel(sprintf('frequency (difference)'));gcaf1(gca,35);
setAxesForPPT(ax1, 35, 2.5);cutBackColor(gca,gcf);
print(gcf,[f_file,'modelGamma_diff_weak_strong.svg'], '-dsvg','-vector','-r600');
%% figure3O
time_all=double(t_wind_start+100)./1000;
freq_all=f;
time_index=find(round(time_all,3)==-0.5):find(round(time_all,3)==-0.3);
freq_index=round(freq_all)>=40 & round(freq_all)<70;
for i=1:size(low_mtfr,1)
os_vtfr(i,1)=squeeze(mean(low_mtfr(i,freq_index,time_index),'all'));
ps_vtfr(i,1)=squeeze(mean(high_mtfr(i,freq_index,time_index),'all'));
end
xx1=[os_vtfr(order_idx1),ps_vtfr(order_idx2)];
ax1=figureX1([24 20]);
x_pos=figureTmp1_OT_SinglePair(xx1,1);
box off;ylabel(sprintf('gamma-band power (log)'));
set(gca,'XTick',x_pos,'XTickLabel',{'weak','strong'}); xlabel('mutual inhibition (λ)');
setFixedPlotArea(ax1, [15 15]);gcaf1(gca,35);setAxesForPPT(ax1, 35, 2.5);
cutBackColor(gca,gcf);print(gcf,[f_file,'modelGamma_tfr_interaction(40-70hz)'], '-dsvg','-vector','-r600');
%% figure3P
int_color = [132 121 181]./255;int_color1=[132 121 181]./255;
xx2=os_vtfr(order_idx1)-ps_vtfr(order_idx2);
lambda_int=-inhibit_low(order_idx1)'+inhibit_high(order_idx2)';
ax1=figureX1([24 20]);
hold on;
scatter(lambda_int,xx2,'MarkerEdgeColor',[1 1 1],'MarkerFaceColor',int_color,'SizeData',305,'marker','square');
pl=errorCIbound(lambda_int,xx2,int_color1,linspace(-0.1,0,100)');pl.LineWidth=3;
hold off;xlim([-0.1 0]);
xlabel('inhibition difference (Δλ)');ylabel(sprintf('gamma-band power difference (log)'));
setFixedPlotArea(ax1, [15 15]);set(gca, 'LineWidth', 3, 'TickDir', 'in'); 
get2CorrInfor(lambda_int,xx2);
setAxesForPPT(ax1, 35, 2.5);gcaf1(gca,35);
print(gcf,[f_file,'corr(model_lambda_int,tfr_13-30hz)'], '-dsvg','-vector','-r600');


%% figure3 Q-T
% MEG
xx_tmp=load(fullfile(dataRoot,'MEG_results','high_powspectrm_congruentpressPutton_allCondition(250hz)102.mat'));
ppm_cg.all=xx_tmp.(string(fields(xx_tmp)));
xx_tmp=load(fullfile(dataRoot,'MEG_results','high_powspectrm_conflictpressPutton_allCondition(250hz)102.mat'));
ppm_cf.all=xx_tmp.(string(fields(xx_tmp)));

% figure3Q
% sensor level: channel 23
target_channel_idx = 23;
target_channel_label = ppm_cf.all.label{target_channel_idx};
cfg_sel = [];
cfg_sel.channel = target_channel_label;
condition1 = ft_selectdata(cfg_sel, ppm_cf.all);
condition2 = ft_selectdata(cfg_sel, ppm_cg.all);
condition_d = condition2;condition_d.powspctrm=log(condition1.powspctrm)-log(condition2.powspctrm);
cfg=[];
cfg.frequency_window = [30 100];
cfg.latency_window   = [-1 0];
cfg.clim=[-4 4];
TFR_figure(condition_d,cfg);cutBackColor(gca,gcf);
print(gcf,[f_file,'channel23_time-frequency_Tvalue(30-120hz)'], '-dsvg','-vector','-r600');
%% figure3R
% channel 23 time course
time_index=find(ppm_cf.all.time==-2.0):find(ppm_cf.all.time==0.5);
x_axis=ppm_cf.all.time(time_index);
freq_index=find(round(ppm_cf.all.freq)==60):find(round(ppm_cf.all.freq)==80);
xx_int=cell(size(ppm_cf.all.powspctrm,2),1);
for i_channel=1:size(ppm_cf.all.powspctrm,2)
    xx_int{i_channel,1}=squeeze(mean(log(ppm_cf.all.powspctrm(:,i_channel,freq_index,time_index))-log(ppm_cg.all.powspctrm(:,i_channel,freq_index,time_index)),[3]));
end
%
int_color = [132 121 181]./255;
ax1=figureX1([24 20]);
i_channel=23;
hold on
errorbound(x_axis,mean(xx_int{i_channel}),std(xx_int{i_channel})./sqrt(size(xx_int{i_channel},1)),'color',int_color,'linewidth',5);
plot_column_significance(x_axis,xx_int{i_channel},'Color',int_color, 'LineWidth', 5,'YOffsetRatio',0.04,'TimeWindow',[-0.8 0]);
xline(0,'k--','LineWidth',2.5);
yline(0,'k--','LineWidth',2.5);
setFixedPlotArea(ax1, [15 15]);set(gca, 'LineWidth', 1.5); 
xlim([-1 0.5]);
hold off;ylabel(sprintf('power(log)\n conflict—congruent'));xlabel('time from choice made');
gcaf1(gca,35);setAxesForPPT(ax1, 35, 2.5);cutBackColor(gca,gcf);
print(gcf,[f_file,'time-course_(60-80hz)'], '-dsvg','-vector','-r600');
%% figure3S
% topoplot Tvalue
load(fullfile(dataRoot,'MEG_results','fre_60_80Hz_topoplot.mat'));
ot_freq=cell(length(ot_freq_conf),1);pl_freq=cell(length(ot_freq_conf),1);all_freq=cell(length(ot_freq_conf),1);
for i=1:length(ot_freq_conf)
ot_freq{i}=ot_freq_conf{i};
ot_freq{i}.powspctrm=log(ot_freq_conf{i}.powspctrm)-log(ot_freq_cong{i}.powspctrm);
pl_freq{i}=pl_freq_conf{i};
pl_freq{i}.powspctrm=log(pl_freq_conf{i}.powspctrm)-log(pl_freq_cong{i}.powspctrm);
all_freq{i}=ot_freq_conf{i};
all_freq{i}.powspctrm=(ot_freq{i}.powspctrm+pl_freq{i}.powspctrm)./2;
end
cfg_grand = [];
cfg_grand.keepindividual = 'yes';
tmp_freq0 = ft_freqgrandaverage(cfg_grand, all_freq{:});
[mask_tmp,~,~,stat_tmp]=ttest(tmp_freq0.powspctrm);
tmp_freq=tmp_freq0;tmp_freq.powspctrm=stat_tmp.tstat';
tmp_freq.sig_mask=mask_tmp';tmp_freq.sig_mask(mask_tmp==1)=1;tmp_freq.sig_mask(mask_tmp==0)=0.18;
%
cfg = [];
cfg.xlim         =[-0.5 -0.3];
cfg.zlim         = [-3 3];
cfg.marker       = 'off';
cfg.layout       = 'neuromag306cmb.lay';
cfg.comment      = 'no';
cfg.colorbar     = 'no';
cfg.highlight  ='on';
cfg.maskparameter = 'sig_mask';
cfg.maskstyle     = 'opacity';
cfg.maskalpha     = tmp_freq.sig_mask;    
cfg.highlightchannel = tmp_freq.label([23]);
cfg.highlightsymbol = '*';
cfg.highlightcolor  = [0 0 0];   % Red
cfg.highlightsize   = 24;
cmap = get_adaptive_cmap(cfg.zlim);
cfg.figure=figureX1([23 20]);
ft_topoplotER(cfg, tmp_freq);
colormap(cmap);
setFixedPlotArea(cfg.figure, [15 15]);
print(gcf,[f_file,'choice_60-80hz_topgraphic——Tvalue(source_use)'], '-dsvg','-vector','-r600');

%% figure3T
% brainNet to generate gamma-band source ∩ oxytocin reciptor gene
% expression map
cfg = struct();
cfg.brainnet_dir = fullfile(codeRoot,'utilities','BrainNet-Viewer');
cfg.surface_file = fullfile(cfg.brainnet_dir,'BrainMesh_ICBM152Left.nv');
cfg.cfg_file = fullfile(dataRoot,'MEG_results','BrainNet_cfg_medial_left.mat');

cfg.mapping_dir = fullfile(dataRoot,'MEG_results');
cfg.output_dir = f_file;

% Render the map assigned to this panel.
cfg.input_pattern = 'decision_gamma_band_source_OXTRconjunction_.nii';

% Preserve saved display thresholds and color limits.
cfg.fix_global_color_scale = false;
cfg.save_tif_only = true;
result = batch_brainnet_nii_to_tif(cfg);

%% figure 3U-X
% load data for figure3U-X
load(fullfile(dataRoot,'MEG_results','decision_gamma_band_sourceTFR.mat'));
[high_all_tfr_diff,high_os_tfr_diff, ...
 high_ps_tfr_diff,high_on_tfr_diff,high_pn_tfr_diff] = ...
    getConstractPower(all_tfr_conf,all_tfr_cong,os_tfr_conf,os_tfr_cong, ...
        ps_tfr_conf,ps_tfr_cong,on_tfr_conf,on_tfr_cong,pn_tfr_conf,pn_tfr_cong);

cfg_ga1 = [];
cfg_ga1.keepindividual = 'yes';
all_tfr_diff1 = ft_freqgrandaverage(cfg_ga1, high_all_tfr_diff{:});
os_tfr_diff1 = ft_freqgrandaverage(cfg_ga1, high_os_tfr_diff{:});
ps_tfr_diff1 = ft_freqgrandaverage(cfg_ga1, high_ps_tfr_diff{:});
on_tfr_diff1 = ft_freqgrandaverage(cfg_ga1, high_on_tfr_diff{:});
pn_tfr_diff1 = ft_freqgrandaverage(cfg_ga1, high_pn_tfr_diff{:});


%% figure3U
% source single TFR: average: gamma
int_so=all_tfr_diff1;
cfg=[];
cfg.frequency_window=[50 70];
cfg.latency_window=[-0.5 -0.3];
cfg.clim=[-3 1.5];
TFR_figure(int_so,cfg);
print(gcf,[f_file,'source_all_time-frequency_Tvalue(30-120hz)'], '-dsvg','-vector','-r600');
%% figure3V
% source single TFR: interaction: gamma
int_so=os_tfr_diff1;
int_so.powspctrm = os_tfr_diff1.powspctrm - ps_tfr_diff1.powspctrm - on_tfr_diff1.powspctrm + pn_tfr_diff1.powspctrm;
cfg=[];
cfg.frequency_window=[50 70];
cfg.latency_window=[-0.5 -0.3];
cfg.clim=[-3 3];
TFR_figure(int_so,cfg);
print(gcf,[f_file,'source_int_time-frequency_Tvalue(30-120hz)_allSigcluster'], '-dsvg','-vector','-r600');

%% figure3U
% 40-70Hz power interaction
time_all=all_tfr_diff1.time;
freq_all=all_tfr_diff1.freq;
time_index=(round(time_all,3)>=-0.5) & (round(time_all,3)<=-0.3);
freq_index=(round(freq_all)>40) & (round(freq_all)<=70);
for i=1:length(high_os_tfr_diff)
os_vtfr(i,1)=squeeze(mean(high_os_tfr_diff{i}.powspctrm(1,freq_index,time_index),'all'));
ps_vtfr(i,1)=squeeze(mean(high_ps_tfr_diff{i}.powspctrm(1,freq_index,time_index),'all'));
on_vtfr(i,1)=squeeze(mean(high_on_tfr_diff{i}.powspctrm(1,freq_index,time_index),'all'));
pn_vtfr(i,1)=squeeze(mean(high_pn_tfr_diff{i}.powspctrm(1,freq_index,time_index),'all'));
end
%
xx1=[os_vtfr,ps_vtfr,on_vtfr,pn_vtfr];
ax1=figureX1([24 20]);
x_pos=figureTmp1_OT(xx1,1);
box off;ylabel(sprintf('power(log)\n conflict—congruent'));
set(gca,'XTick',x_pos,'XTickLabel',{'OT','PL','OT','PL'}); gcaf1(gca,35);
setFixedPlotArea(ax1, [15 15]);setAxesForPPT(ax1, 35, 2.5);
cutBackColor(gca,gcf);print(gcf,[f_file,'tfr_interaction(40-70hz)'], '-dsvg','-vector','-r600');
%% figure3X
% corr (λ interaction , power interaction)
int_color = [132 121 181]./255;int_color1=[132 121 181]./255;
xx2=os_vtfr-on_vtfr-ps_vtfr+pn_vtfr;
lambda_int=[ots_xx1(:,7)-otn_xx1(:,7)-pls_xx1(:,7)+pln_xx1(:,7)];
ax1=figureX1([24 20]);
hold on;
scatter(lambda_int,xx2,'MarkerEdgeColor',[1 1 1],'MarkerFaceColor',int_color,'SizeData',305,'marker','square');
pl=errorCIbound(lambda_int,xx2,int_color1,[-0.2:0.01:0.1]');pl.LineWidth=3;
hold off;
xlabel('interaction effect (λ)');ylabel(sprintf('gamma-band power (log)\n(interaction)'));gcaf1(gca,35);
setFixedPlotArea(ax1, [15 15],[1.2 0.6]);xlim tight;set(gca, 'LineWidth', 3, 'TickDir', 'in'); 
setAxesForPPT(ax1, 35, 2.5);
get2CorrInfor(lambda_int,xx2);
cutBackColor(gca,gcf);print(gcf,[f_file,'corr(lambda_int,tfr_40-70hz)'], '-dsvg','-vector','-r600');














%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [high_all_tfr_differ,high_os_tfr_differ, ...
    high_ps_tfr_differ,high_on_tfr_differ,high_pn_tfr_differ] = ...
    getConstractPower(all_conf,all_cong,os_conf,os_cong, ...
    ps_conf,ps_cong,on_conf,on_cong,pn_conf,pn_cong)
% Outputs: subject-by-1 cells containing log(conflict)-log(congruent) TFR.
% Inputs must have matching subject order and channel/frequency/time grids.
high_all_tfr_differ = log_contrast_local(all_conf,all_cong);
high_os_tfr_differ  = log_contrast_local(os_conf,os_cong);
high_ps_tfr_differ  = log_contrast_local(ps_conf,ps_cong);
high_on_tfr_differ  = log_contrast_local(on_conf,on_cong);
high_pn_tfr_differ  = log_contrast_local(pn_conf,pn_cong);
end

function result = log_contrast_local(conf,cong)
result = cell(numel(conf),1);
for i = 1:numel(conf)
    result{i} = conf{i};
    result{i}.powspctrm = log(conf{i}.powspctrm) ...
                       - log(cong{i}.powspctrm);
    if isfield(result{i},'cfg')
        result{i} = rmfield(result{i},'cfg');
    end
end
end

function TFR_figure(int_so,cfg_test)
int_so0=int_so;int_so0.powspctrm=zeros(size(int_so.powspctrm,1),size(int_so.powspctrm,2),size(int_so.powspctrm,3),size(int_so.powspctrm,4));
% F_pointwise=getF_test_noCluster(int_so,int_so0);
F_pointwise=getF_test_cluster_size(int_so,int_so0,cfg_test.frequency_window,cfg_test.latency_window);
xx_tmp=getT_value(int_so);
int_t_value=int_so;int_t_value.powspctrm=xx_tmp;int_t_value.dimord='chan_freq_time';
int_t_value = attach_cluster_mask_to_tfr(int_t_value,F_pointwise,'cluster_mask');
int_smooth=getSmoothTFR(int_t_value);
channel_label = int_smooth.label{1};
channel_idx = find(strcmp(int_smooth.label,channel_label),1);
fig = figureX1([23,20]);
cfg = [];
cfg.figure        = fig;
cfg.parameter     = 'powspctrm';
cfg.channel       = channel_label;
cfg.title         = ' ';
cfg.colorbar      = 'no';
cfg.feedback      = 'none';
cfg.xlim = [-1 0.5];
cfg.zlim = cfg_test.clim;%[-3 3];
cfg.maskparameter = 'cluster_mask';
cfg.maskstyle = 'opacity';
% Opacity of nonsignificant regions
cfg.maskalpha = 1;
cmap = get_adaptive_cmap(cfg.zlim);
ft_singleplotTFR(cfg,int_smooth);
ax = gca;colormap(ax,cmap);set(ax,'LineWidth',2);
hold(ax,'on');
sig_line(ax,int_smooth,channel_idx);
xline(ax,0,'k--','LineWidth',2.5);
xlabel(ax,'time from choice made');
ylabel(ax,'frequency (Hz)');
xlim(ax,[-1 0.5]);
setFixedPlotArea(fig,[15 15]);
gcaf1(ax,35);hold(ax,'off');
setAxesForPPT(fig, 35, 2.5);
end


function F_cluster = getF_test_cluster_size(condition1, condition2, frequency_window, latency_window)

%% Check data
assert(isequal(condition1.label, condition2.label), ...
    '两个条件的通道标签不一致。');

assert(isequal(condition1.freq, condition2.freq), ...
    '两个条件的频率轴不一致。');

assert(isequal(condition1.time, condition2.time), ...
    '两个条件的时间轴不一致。');

n_subject = size(condition1.powspctrm,1);

assert(size(condition2.powspctrm,1) == n_subject, ...
    '两个条件的被试数不一致。');

%% Within-participant design matrix
design = zeros(2,2*n_subject);

design(1,:) = [ones(1,n_subject), ...
               2*ones(1,n_subject)];

design(2,:) = [1:n_subject, ...
               1:n_subject];

%% Cluster permutation F test
cfg_stat = [];

cfg_stat.parameter = 'powspctrm';
cfg_stat.channel   = 'all';
cfg_stat.frequency = frequency_window;
cfg_stat.latency   = latency_window;

cfg_stat.method    = 'montecarlo';
cfg_stat.statistic = 'ft_statfun_depsamplesFunivariate';

cfg_stat.correctm = 'cluster';

% Cluster Forming threshold, Not the final significance threshold
cfg_stat.clusteralpha = 0.05;

% Key point: Use cluster Magnitude as the cluster statistic
cfg_stat.clusterstatistic = 'maxsize';%'maxsum';%%'maxsum';'max''wcm';

% No spatial adjacency for a single channel
cfg_stat.neighbours = [];
cfg_stat.minnbchan  = 0;

% F values are positive only
cfg_stat.tail        = 1;
cfg_stat.clustertail = 1;

% Final cluster significance threshold
cfg_stat.alpha = 0.05;

cfg_stat.numrandomization = 1000;

cfg_stat.design = design;
cfg_stat.ivar   = 1;
cfg_stat.uvar   = 2;

tic;
F_cluster = ft_freqstatistics( ...
    cfg_stat, condition1, condition2);
elapsed_time = toc;

%% Organize outputs
F_cluster.F_all = F_cluster.stat;

F_cluster.F_sig = F_cluster.stat;
F_cluster.F_sig(~F_cluster.mask) = 0;

fprintf('\n===================================================\n');
fprintf('Repeated-measures cluster F test\n');
fprintf('Channel           : %s\n', condition1.label{1});
fprintf('Frequency         : [%g %g] Hz\n', frequency_window);
fprintf('Latency           : [%g %g] s\n', latency_window);
fprintf('Subjects          : %d\n', n_subject);
fprintf('Cluster statistic : maxsize\n');
fprintf('Randomizations    : %d\n', cfg_stat.numrandomization);
fprintf('Significant bins  : %d\n', nnz(F_cluster.mask));
fprintf('Elapsed time      : %.2f s\n', elapsed_time);

if isfield(F_cluster,'posclusters') && ...
        ~isempty(F_cluster.posclusters)

    for i_cluster = 1:numel(F_cluster.posclusters)

        cluster_size = nnz( ...
            F_cluster.posclusterslabelmat == i_cluster);

        fprintf('Cluster %d: size=%d, p=%.6f\n', ...
            i_cluster, ...
            cluster_size, ...
            F_cluster.posclusters(i_cluster).prob);
    end
else
    fprintf('No positive cluster was found.\n');
end

fprintf('===================================================\n');
end

function sig_line(ax,int_smooth,channel_idx)
mask_2d = squeeze(int_smooth.cluster_mask(channel_idx,:,:));
if any(mask_2d(:))
    contour(ax, ...
        int_smooth.time, ...
        int_smooth.freq, ...
        double(mask_2d), ...
        [0.5 0.5], ...
        'LineColor',[0.7 0.7 0.7], ...
        'LineWidth',2.5, ...
        'linestyle','--');
end
end

function xx_tmp=getT_value(int_tfr_diff1)
% subj_chan_freq_time
xx_tmp=nan(length(int_tfr_diff1.label),length(int_tfr_diff1.freq),length(int_tfr_diff1.time));
for i_channel=1:length(int_tfr_diff1.label)
    for i=1:length(int_tfr_diff1.freq)
        for j=1:length(int_tfr_diff1.time)
            xx1=int_tfr_diff1.powspctrm(:,i_channel,i,j);
            if sum(isnan(xx1))==0
                [~,~,~,stat]=ttest(xx1);
                xx_tmp(i_channel,i,j)=stat.tstat;
            else
                xx_tmp(i_channel,i,j)=nan;
            end
        end
    end
end
end


function tfr_plot = attach_cluster_mask_to_tfr(tfr_plot, stat, mask_field)
% ATTACH_CLUSTER_MASK_TO_TFR Align a FieldTrip cluster mask to a TFR grid.
%
% tfr_plot must contain chan_freq_time data. stat is the output of
% ft_freqstatistics and must contain stat.mask, stat.label, stat.freq and
% stat.time. Frequencies/times tested in stat are inserted into the matching
% locations of the potentially larger plotting TFR; all untested locations
% remain false.

if nargin < 3 || isempty(mask_field)
    mask_field = 'cluster_mask';
end

assert(isfield(tfr_plot, 'powspctrm') && ...
       strcmp(tfr_plot.dimord, 'chan_freq_time'), ...
    'tfr_plot must have dimord = chan_freq_time.');
assert(isfield(stat, 'mask') && isfield(stat, 'label') && ...
       isfield(stat, 'freq') && isfield(stat, 'time'), ...
    'stat must contain mask, label, freq and time.');

stat_mask = logical(stat.mask);
expected_size = [numel(stat.label), numel(stat.freq), numel(stat.time)];
stat_mask = reshape(stat_mask, expected_size);

% Map the statistics channels onto the plotting channels.
[channel_found, channel_index] = ismember(stat.label, tfr_plot.label);
assert(all(channel_found), ...
    'At least one statistics channel is absent from tfr_plot.label.');

% Match numerical axes using a tolerance rather than exact floating-point
% equality.
freq_index = match_numeric_axis(stat.freq, tfr_plot.freq, 'frequency');
time_index = match_numeric_axis(stat.time, tfr_plot.time, 'time');

plot_mask = false(size(tfr_plot.powspctrm));
plot_mask(channel_index, freq_index, time_index) = stat_mask;
tfr_plot.(mask_field) = plot_mask;
end

function target_index = match_numeric_axis(source_axis, target_axis, axis_name)
source_axis = double(source_axis(:));
target_axis = double(target_axis(:));
target_index = nan(size(source_axis));

if numel(target_axis) > 1
    tolerance = max(1e-10, min(abs(diff(target_axis))) * 1e-5);
else
    tolerance = 1e-10;
end

for k = 1:numel(source_axis)
    [distance, idx] = min(abs(target_axis - source_axis(k)));
    assert(distance <= tolerance, ...
        'Cannot match statistics %s value %.12g to the plotting axis.', ...
        axis_name, source_axis(k));
    target_index(k) = idx;
end

target_index = target_index(:).';
end


% ================= 1. Smoothly interpolate original time-frequency data =================
function TFR_smooth = getSmoothTFR(TFR_tmp1)
% =============================================================
% For TFR power maps and statistical mask use different interpolation types: 
%
% powspctrm   : spline Interpolation, For smooth display
% cluster_mask: nearest Interpolation, Preserve binary statistical results
%
% Note: Interpolated mask For plotting only, Do not reuse for statistical testing.
% =============================================================

assert(strcmp(TFR_tmp1.dimord,'chan_freq_time'), ...
    '输入数据必须为chan_freq_time。');

time_orig = double(TFR_tmp1.time(:)');
freq_orig = double(TFR_tmp1.freq(:)');

n_time_orig = numel(time_orig);
n_freq_orig = numel(freq_orig);
n_channel   = numel(TFR_tmp1.label);

assert(isequal(size(TFR_tmp1.powspctrm), ...
    [n_channel,n_freq_orig,n_time_orig]), ...
    'powspctrm维度与label/freq/time不一致。');

% -------------------------------------------------------------
% Interpolation factor
% -------------------------------------------------------------
interp_factor = 2;

% Use (n-1)*factor+1, Ensure original samples remain in the new grid
n_time_fine = ...
    (n_time_orig-1)*interp_factor + 1;

n_freq_fine = ...
    (n_freq_orig-1)*interp_factor + 1;

time_fine = linspace( ...
    time_orig(1),time_orig(end),n_time_fine);

freq_fine = linspace( ...
    freq_orig(1),freq_orig(end),n_freq_fine);

[T_orig,F_orig] = meshgrid(time_orig,freq_orig);
[T_fine,F_fine] = meshgrid(time_fine,freq_fine);

% -------------------------------------------------------------
% Initialize output structure
% -------------------------------------------------------------
TFR_smooth = TFR_tmp1;

TFR_smooth.time = time_fine;
TFR_smooth.freq = freq_fine;

TFR_smooth.powspctrm = nan( ...
    n_channel,n_freq_fine,n_time_fine);

% -------------------------------------------------------------
% Power maps use spline Interpolation
% -------------------------------------------------------------
for i_chan = 1:n_channel

    pow_orig = squeeze( ...
        TFR_tmp1.powspctrm(i_chan,:,:));

    % If this channel exists NaN, Use instead linear, Prevent spline expansion of NaN region
    if any(~isfinite(pow_orig(:)))
        interpolation_method = 'linear';
    else
        interpolation_method = 'spline';
    end

    pow_fine = interp2( ...
        T_orig,F_orig,pow_orig, ...
        T_fine,F_fine, ...
        interpolation_method);

    TFR_smooth.powspctrm(i_chan,:,:) = pow_fine;
end

% -------------------------------------------------------------
% significance mask Use nearest Interpolation
% -------------------------------------------------------------
if isfield(TFR_tmp1,'cluster_mask')

    mask_orig_all = logical(TFR_tmp1.cluster_mask);

    assert(isequal(size(mask_orig_all), ...
        size(TFR_tmp1.powspctrm)), ...
        ['cluster_mask必须与原始powspctrm维度一致。', ...
         '请先使用attach_cluster_mask_to_tfr进行对齐。']);

    TFR_smooth.cluster_mask = false( ...
        n_channel,n_freq_fine,n_time_fine);

    for i_chan = 1:n_channel

        mask_orig = squeeze( ...
            mask_orig_all(i_chan,:,:));

        % Nearest-neighbor interpolation, Does not generate statistical values between 0 and 1 these values
        mask_fine = interp2( ...
            T_orig,F_orig,double(mask_orig), ...
            T_fine,F_fine, ...
            'nearest');

        TFR_smooth.cluster_mask(i_chan,:,:) = ...
            mask_fine >= 0.5;
    end
end

% Record that this structure is for smooth display only
TFR_smooth.display_interpolation = struct();
TFR_smooth.display_interpolation.factor = interp_factor;
TFR_smooth.display_interpolation.power_method = ...
    'spline_or_linear_if_nan';
TFR_smooth.display_interpolation.mask_method = ...
    'nearest';
TFR_smooth.display_interpolation.statistical_use = ...
    'display_only';

end

%% function 
function cmap = get_adaptive_cmap(zlim_vals, n_colors, use_input_scale)
% GET_ADAPTIVE_CMAP Generate adaptive coolwarm style colormap.
%
% Input parameters:
%   zlim_vals       : 1x2 vector, color/Z axis range, For example [-3.5, 5.2]
%   n_colors        : (Optional) colormap Order, Default 256
%   use_input_scale : (Optional) colormap Mapping mode: 
%                     Empty value or 0(Default)——Preserve original logic, When spanning 0 set 0
%                     strictly to neutral gray-white; 
%                     1 —— Complete coolwarm Color scale covers the input range, Gray-white lies at
%                     zlim_vals the midpoint, Do not force numeric 0 to gray-white.
%
% Output parameters:
%   cmap            : n_colors x 3 's RGB matrix

    if nargin < 2 || isempty(n_colors)
        n_colors = 256;
    end
    if nargin < 3 || isempty(use_input_scale)
        use_input_scale = 0;
    end

    if ~isnumeric(zlim_vals) || numel(zlim_vals) ~= 2 || ...
            any(~isfinite(zlim_vals))
        error('输入错误：zlim_vals 必须是包含两个有限数值的向量。');
    end
    if ~isnumeric(n_colors) || ~isscalar(n_colors) || ...
            ~isfinite(n_colors) || n_colors < 2 || n_colors ~= fix(n_colors)
        error('输入错误：n_colors 必须是大于等于 2 的整数。');
    end
    if ~(isnumeric(use_input_scale) || islogical(use_input_scale)) || ...
            ~isscalar(use_input_scale) || ~ismember(use_input_scale, [0, 1])
        error('输入错误：use_input_scale 必须为空、0 或 1。');
    end

    zmin = zlim_vals(1);
    zmax = zlim_vals(2);

    if zmin >= zmax
        error('输入错误：zlim_vals(1) 必须严格小于 zlim_vals(2)。');
    end

    % 7 level custom version coolwarm Color anchor points(Taken from Python matplotlib style).
    c_dark_blue  = [ 59,  76, 192] / 255;
    c_mid_blue   = [103, 136, 238] / 255;
    c_light_blue = [154, 187, 255] / 255;
    c_center     = [221, 221, 221] / 255;
    c_light_red  = [245, 158, 114] / 255;
    c_mid_red    = [217,  88,  71] / 255;
    c_dark_red   = [180,   4,  38] / 255;

    if use_input_scale == 1
        % Spread the full color scale uniformly over the input range; Neutral color corresponds to the input midpoint.
        vals = linspace(zmin, zmax, 7);
        colors = [c_dark_blue; c_mid_blue; c_light_blue; c_center; ...
                  c_light_red; c_mid_red; c_dark_red];
    elseif zmin < 0 && zmax > 0
        % Spanning zero: Allocate colors according to positive and negative ranges, Ensure 0 Corresponds to neutral gray-white.
        vals = [zmin, zmin * 0.66, zmin * 0.33, 0, ...
                zmax * 0.33, zmax * 0.66, zmax];
        colors = [c_dark_blue; c_mid_blue; c_light_blue; c_center; ...
                  c_light_red; c_mid_red; c_dark_red];
    elseif zmin >= 0
        % All positive or including 0: Gray-white -> Red.
        vals = linspace(zmin, zmax, 4);
        colors = [c_center; c_light_red; c_mid_red; c_dark_red];
    else
        % All negative: Blue -> Gray-white.
        vals = linspace(zmin, zmax, 4);
        colors = [c_dark_blue; c_mid_blue; c_light_blue; c_center];
    end

    vals_norm = (vals - zmin) / (zmax - zmin);
    query_points = linspace(0, 1, n_colors);
    cmap = interp1(vals_norm, colors, query_points, 'pchip');
    cmap = max(min(cmap, 1), 0);
end


%%

function h = plot_gradient_trajectory_sem( ...
        subjectX,subjectY,lightColor,darkColor,options)
%PLOT_GRADIENT_TRAJECTORY_SEM Plot a 2-D mean path with a paired SEM band.
% subjectX and subjectY are subject-by-step matrices. At each step, paired
% subject deviations are projected onto the local normal of the mean path.
% The ribbon is mean path +/- semMultiplier times the normal-direction SEM.
%
% Example:
%   h = plot_gradient_trajectory_sem(staySelf,stayOther, ...
%       lightColor,darkColor,lineWidth=4,semAlpha=0.20, ...
%       startMarker='square',endMarker='o', ...
%       startMarkerSize=120,endMarkerSize=120);

    arguments
        subjectX (:,:) double
        subjectY (:,:) double
        lightColor (1,3) double
        darkColor (1,3) double
        options.semMultiplier (1,1) double = 1
        options.semAlpha (1,1) double = 0.18
        options.semColorStrength (1,1) double = 0.45
        options.minimumSubjects (1,1) double = 3
        options.lineWidth (1,1) double = 5
        options.showStartMarker (1,1) logical = true
        options.showEndMarker (1,1) logical = true
        options.startMarker (1,:) char = 'square'
        options.endMarker (1,:) char = 'o'
        options.startMarkerSize (1,1) double = 180
        options.endMarkerSize (1,1) double = 150
        options.markerLineWidth (1,1) double = 1.5
        options.markerFaceColor (1,3) double = [1 1 1]
        options.tangentSmoothWindow (1,1) double = 15
    end

    assert(isequal(size(subjectX),size(subjectY)), ...
        'subjectX and subjectY must have identical dimensions.');
    assert(size(subjectX,1) >= 2 && size(subjectX,2) >= 2, ...
        'Input must contain at least two subjects and two steps.');
    validate_color_local(lightColor,'lightColor');
    validate_color_local(darkColor,'darkColor');
    validate_color_local(options.markerFaceColor,'markerFaceColor');
    assert(options.semMultiplier > 0 && isfinite(options.semMultiplier), ...
        'semMultiplier must be positive and finite.');
    assert(options.semAlpha >= 0 && options.semAlpha <= 1, ...
        'semAlpha must lie between zero and one.');
    assert(options.semColorStrength >= 0 && ...
        options.semColorStrength <= 1, ...
        'semColorStrength must lie between zero and one.');
    assert(options.minimumSubjects >= 2 && ...
        options.minimumSubjects == round(options.minimumSubjects), ...
        'minimumSubjects must be an integer of at least two.');
    assert(options.lineWidth > 0 && isfinite(options.lineWidth), ...
        'lineWidth must be positive and finite.');
    assert(options.startMarkerSize > 0 && options.endMarkerSize > 0, ...
        'Marker sizes must be positive.');
    assert(options.markerLineWidth >= 0 && ...
        isfinite(options.markerLineWidth), ...
        'markerLineWidth must be nonnegative and finite.');

    subjectX = double(subjectX);
    subjectY = double(subjectY);
    pairValid = isfinite(subjectX) & isfinite(subjectY);
    subjectX(~pairValid) = nan;
    subjectY(~pairValid) = nan;
    nValid = sum(pairValid,1);
    supported = nValid >= options.minimumSubjects;
    meanX = mean(subjectX,1,'omitnan');
    meanY = mean(subjectY,1,'omitnan');
    meanX(~supported) = nan;
    meanY(~supported) = nan;
    nPoint = numel(meanX);

    [normalX,normalY] = trajectory_normals_local( ...
        meanX,meanY,supported,options.tangentSmoothWindow);
    normalDeviation = (subjectX-meanX).*normalX + ...
        (subjectY-meanY).*normalY;
    semNormal = std(normalDeviation,0,1,'omitnan')./sqrt(nValid);
    semNormal(~supported) = nan;
    halfWidth = options.semMultiplier*semNormal;

    lowerX = meanX-halfWidth.*normalX;
    lowerY = meanY-halfWidth.*normalY;
    upperX = meanX+halfWidth.*normalX;
    upperY = meanY+halfWidth.*normalY;
    semBand = union_local_quadrilaterals( ...
        lowerX,lowerY,upperX,upperY);

    progress = linspace(0,1,nPoint)';
    lineColors = (1-progress).*lightColor + progress.*darkColor;
    middleColor = 0.5*(lightColor+darkColor);
    semColor = (1-options.semColorStrength).*ones(1,3) + ...
        options.semColorStrength.*middleColor;

    ax = gca;
    hold(ax,'on');
    if isempty(semBand.Vertices)
        h.sem = gobjects(0);
    else
        h.sem = plot(ax,semBand,'FaceColor',semColor, ...
            'FaceAlpha',options.semAlpha,'EdgeColor','none', ...
            'HandleVisibility','off');
    end

    surfaceColors = reshape(lineColors,1,nPoint,3);
    surfaceColors = repmat(surfaceColors,2,1,1);
    h.line = surface(ax,[meanX;meanX],[meanY;meanY], ...
        zeros(2,nPoint),surfaceColors,'FaceColor','none', ...
        'EdgeColor','interp','LineWidth',options.lineWidth, ...
        'HandleVisibility','off');

    finiteMean = find(isfinite(meanX) & isfinite(meanY));
    h.start = gobjects(0);
    h.end = gobjects(0);
    if ~isempty(finiteMean) && options.showStartMarker
        first = finiteMean(1);
        h.start = scatter(ax,meanX(first),meanY(first), ...
            options.startMarkerSize,'Marker',options.startMarker, ...
            'MarkerFaceColor',options.markerFaceColor, ...
            'MarkerEdgeColor',lightColor, ...
            'LineWidth',options.markerLineWidth, ...
            'HandleVisibility','off');
    end
    if ~isempty(finiteMean) && options.showEndMarker
        last = finiteMean(end);
        h.end = scatter(ax,meanX(last),meanY(last), ...
            options.endMarkerSize,'Marker',options.endMarker, ...
            'MarkerFaceColor',darkColor, ...%options.markerFaceColor
            'MarkerEdgeColor',darkColor, ...
            'LineWidth',options.markerLineWidth, ...
            'HandleVisibility','off');
    end

    h.mean.self = meanX;
    h.mean.other = meanY;
    h.semNormal = semNormal;
    h.semMultiplier = options.semMultiplier;
    h.support.nSubjects = nValid;
    h.support.validStep = supported;
    h.normal.self = normalX;
    h.normal.other = normalY;
    h.options = options;
    axis(ax,'equal');
    box(ax,'off');
end

function [normalX,normalY] = trajectory_normals_local( ...
        meanX,meanY,supported,requestedWindow)
    smoothX = meanX;
    smoothY = meanY;
    valid = isfinite(meanX) & isfinite(meanY) & supported;
    if sum(valid) >= 2
        index = 1:numel(meanX);
        smoothX(~valid) = interp1(index(valid),meanX(valid), ...
            index(~valid),'linear','extrap');
        smoothY(~valid) = interp1(index(valid),meanY(valid), ...
            index(~valid),'linear','extrap');
    end
    smoothWindow = min(round(requestedWindow),numel(meanX));
    if mod(smoothWindow,2) == 0
        smoothWindow = smoothWindow-1;
    end
    if smoothWindow >= 5
        smoothX = smoothdata(smoothX,'sgolay',smoothWindow);
        smoothY = smoothdata(smoothY,'sgolay',smoothWindow);
    end
    dx = gradient(smoothX);
    dy = gradient(smoothY);
    tangentLength = hypot(dx,dy);
    tangentLength(tangentLength < 1e-8) = 1;
    normalX = -dy./tangentLength;
    normalY = dx./tangentLength;
    for point = 2:numel(meanX)
        if dot([normalX(point-1),normalY(point-1)], ...
                [normalX(point),normalY(point)]) < 0
            normalX(point) = -normalX(point);
            normalY(point) = -normalY(point);
        end
    end
    normalX(~valid) = nan;
    normalY(~valid) = nan;
end

function band = union_local_quadrilaterals( ...
        lowerX,lowerY,upperX,upperY)
    validPoint = isfinite(lowerX) & isfinite(lowerY) & ...
        isfinite(upperX) & isfinite(upperY);
    validFace = find(validPoint(1:end-1) & validPoint(2:end));
    band = polyshape();
    for faceIndex = 1:numel(validFace)
        point = validFace(faceIndex);
        x = [lowerX(point),lowerX(point+1), ...
            upperX(point+1),upperX(point)];
        y = [lowerY(point),lowerY(point+1), ...
            upperY(point+1),upperY(point)];
        vertices = unique([x(:),y(:)],'rows','stable');
        if size(vertices,1) < 3
            continue;
        end
        order = convhull(vertices(:,1),vertices(:,2));
        localBand = polyshape(vertices(order,1),vertices(order,2), ...
            'Simplify',true);
        if area(localBand) > 0
            band = union(band,localBand);
        end
    end
end

function validate_color_local(color,name)
    assert(all(isfinite(color)) && all(color >= 0) && all(color <= 1), ...
        '%s must be an RGB vector between zero and one.',name);
end


%%


function out = bootstrap_trajectory_sem( ...
    X,Y,time_idx,varargin)
% X/Y: subject × complete_time
%
% Center trajectory corresponds exactly to: 
% mean across subjects
% → first smooth
% → select time window
% → abs
% → second smooth

p = inputParser;

addParameter(p,'NBoot',5000);
addParameter(p,'FirstSmoothMethod','movmean');
addParameter(p,'FirstSmoothWindow',20);
addParameter(p,'SecondSmoothMethod','gaussian');
addParameter(p,'SecondSmoothWindow',20);
addParameter(p,'Seed',2026);

parse(p,varargin{:});
opt = p.Results;

X = double(X);
Y = double(Y);

assert(isequal(size(X),size(Y)), ...
    'X and Y sizes differ.');

assert(size(X,2)==numel(time_idx), ...
    'time_idx length differs from data time dimension.');

n_sub = size(X,1);
n_selected = nnz(time_idx);

%% Original sample center trajectory
mean_x_signed = mean(X,1,'omitnan');
mean_y_signed = mean(Y,1,'omitnan');

mean_x_signed = smoothdata( ...
    mean_x_signed,2, ...
    opt.FirstSmoothMethod, ...
    opt.FirstSmoothWindow, ...
    'omitnan');

mean_y_signed = smoothdata( ...
    mean_y_signed,2, ...
    opt.FirstSmoothMethod, ...
    opt.FirstSmoothWindow, ...
    'omitnan');

% Identical to the original plotting code: Take absolute values within the selected window, then smooth
mean_x = smoothdata( ...
    abs(mean_x_signed(time_idx)),2, ...
    opt.SecondSmoothMethod, ...
    opt.SecondSmoothWindow, ...
    'omitnan');

mean_y = smoothdata( ...
    abs(mean_y_signed(time_idx)),2, ...
    opt.SecondSmoothMethod, ...
    opt.SecondSmoothWindow, ...
    'omitnan');

%% Paired participants bootstrap
boot_x = nan(opt.NBoot,n_selected,'single');
boot_y = nan(opt.NBoot,n_selected,'single');

rng(opt.Seed,'twister');

for i_boot = 1:opt.NBoot

    % self and other Use exactly the same participant indices
    sub_idx = randi(n_sub,n_sub,1);

    bx = mean(X(sub_idx,:),1,'omitnan');
    by = mean(Y(sub_idx,:),1,'omitnan');

    bx = smoothdata( ...
        bx,2,opt.FirstSmoothMethod, ...
        opt.FirstSmoothWindow,'omitnan');

    by = smoothdata( ...
        by,2,opt.FirstSmoothMethod, ...
        opt.FirstSmoothWindow,'omitnan');

    bx = smoothdata( ...
        abs(bx(time_idx)),2, ...
        opt.SecondSmoothMethod, ...
        opt.SecondSmoothWindow,'omitnan');

    by = smoothdata( ...
        abs(by(time_idx)),2, ...
        opt.SecondSmoothMethod, ...
        opt.SecondSmoothWindow,'omitnan');

    boot_x(i_boot,:) = single(bx);
    boot_y(i_boot,:) = single(by);
end

%% x and y Each one's bootstrap SEM
sem_x = std(boot_x,0,1,'omitnan');
sem_y = std(boot_y,0,1,'omitnan');

%% Local normal to the mean trajectory
% Used only to determine normal directions, Does not change the actual mean trajectory
geometry_n = min(5,n_selected);

geometry_x = smoothdata( ...
    mean_x,2,'movmean',geometry_n,'omitnan');

geometry_y = smoothdata( ...
    mean_y,2,'movmean',geometry_n,'omitnan');

dx = gradient(geometry_x);
dy = gradient(geometry_y);

denom = hypot(dx,dy);

normal_x = nan(1,n_selected);
normal_y = nan(1,n_selected);

valid_normal = isfinite(denom) & denom>eps;

normal_x(valid_normal) = ...
    -dy(valid_normal)./denom(valid_normal);

normal_y(valid_normal) = ...
     dx(valid_normal)./denom(valid_normal);

%% Fill invalid normals
for k = 1:n_selected

    if ~isfinite(normal_x(k)) || ...
            ~isfinite(normal_y(k))

        if k>1
            normal_x(k) = normal_x(k-1);
            normal_y(k) = normal_y(k-1);
        else
            normal_x(k) = 0;
            normal_y(k) = 1;
        end
    end
end

%% Prevent adjacent normals from flipping
for k = 2:n_selected

    if normal_x(k)*normal_x(k-1) + ...
            normal_y(k)*normal_y(k-1)<0

        normal_x(k) = -normal_x(k);
        normal_y(k) = -normal_y(k);
    end
end

%% Bootstrap Project deviations onto the normal
boot_normal = nan(opt.NBoot,n_selected,'single');

for k = 1:n_selected

    boot_normal(:,k) = ...
        (boot_x(:,k)-mean_x(k))*normal_x(k) + ...
        (boot_y(:,k)-mean_y(k))*normal_y(k);
end

% Bootstrap SEM, No further division by sqrt(n_sub)
sem_normal = std(boot_normal,0,1,'omitnan');

%% SEM boundary
upper_x = mean_x + sem_normal.*normal_x;
upper_y = mean_y + sem_normal.*normal_y;

lower_x = mean_x - sem_normal.*normal_x;
lower_y = mean_y - sem_normal.*normal_y;

%% Output
out = struct();

out.mean_x = mean_x;
out.mean_y = mean_y;

out.sem_x = double(sem_x);
out.sem_y = double(sem_y);
out.sem_normal = double(sem_normal);

out.normal_x = normal_x;
out.normal_y = normal_y;

out.upper_x = upper_x;
out.upper_y = upper_y;
out.lower_x = lower_x;
out.lower_y = lower_y;

out.boot_x = boot_x;
out.boot_y = boot_y;

%% SEM Diagnostics relative to the mean trajectory
radius = hypot(mean_x,mean_y);

valid_ratio = ...
    isfinite(radius) & radius>1e-6 & ...
    isfinite(sem_normal);

out.sem_to_radius = nan(size(radius));
out.sem_to_radius(valid_ratio) = ...
    sem_normal(valid_ratio)./radius(valid_ratio);

out.summary = struct();

out.summary.median_sem = ...
    median(sem_normal,'omitnan');

out.summary.max_sem = ...
    max(sem_normal,[],'omitnan');

out.summary.median_radius = ...
    median(radius,'omitnan');

out.summary.median_sem_to_radius = ...
    median(out.sem_to_radius,'omitnan');

out.summary.max_sem_to_radius = ...
    max(out.sem_to_radius,[],'omitnan');

end


function h = draw_trajectory_sem_union(ax,upper_x,upper_y,lower_x,lower_y,sem_color,face_alpha,clip_at_zero)
% DRAW_TRAJECTORY_SEM_UNION Draw a curved 2-D SEM ribbon without alpha buildup.
%
% Each adjacent pair of upper/lower boundary samples defines one local
% quadrilateral. Their geometric union is calculated before plotting, so a
% folded trajectory neither bridges distant parts of the path nor darkens
% where its SEM ribbon overlaps itself.

if nargin<8 || isempty(clip_at_zero)
    clip_at_zero = false;
end

validateattributes(ax,{'matlab.graphics.axis.Axes'},{'scalar'},mfilename,'ax');
upper_x = double(upper_x(:)');
upper_y = double(upper_y(:)');
lower_x = double(lower_x(:)');
lower_y = double(lower_y(:)');
assert(isequal(size(upper_x),size(upper_y),size(lower_x),size(lower_y)), ...
    'All SEM-boundary vectors must have the same size.');
assert(numel(upper_x)>=2,'At least two time samples are required.');
validateattributes(sem_color,{'numeric'}, ...
    {'real','finite','vector','numel',3,'>=',0,'<=',1},mfilename,'sem_color');
validateattributes(face_alpha,{'numeric'}, ...
    {'real','finite','scalar','>=',0,'<=',1},mfilename,'face_alpha');

ribbon_shape = polyshape();
n_time = numel(upper_x);
for k = 1:n_time-1
    vx = [upper_x(k),upper_x(k+1),lower_x(k+1),lower_x(k)];
    vy = [upper_y(k),upper_y(k+1),lower_y(k+1),lower_y(k)];
    if ~all(isfinite([vx,vy]))
        continue;
    end
    segment_shape = polyshape(vx,vy,'Simplify',true);
    if isempty(segment_shape.Vertices) || area(segment_shape)<=eps
        continue;
    end
    if isempty(ribbon_shape.Vertices)
        ribbon_shape = segment_shape;
    else
        ribbon_shape = union(ribbon_shape,segment_shape);
    end
end

if isempty(ribbon_shape.Vertices)
    h = gobjects(0);
    warning('No finite, non-degenerate SEM ribbon could be constructed.');
    return;
end

if clip_at_zero
    vertices = ribbon_shape.Vertices;
    xmax = max(vertices(:,1),[],'omitnan');
    ymax = max(vertices(:,2),[],'omitnan');
    margin = max([xmax,ymax,1])*1e-9;
    clip_shape = polyshape( ...
        [0,xmax+margin,xmax+margin,0], ...
        [0,0,ymax+margin,ymax+margin]);
    ribbon_shape = intersect(ribbon_shape,clip_shape);
end

h = plot(ax,ribbon_shape, ...
    'FaceColor',double(sem_color(:)'), ...
    'FaceAlpha',face_alpha, ...
    'EdgeColor','none', ...
    'HandleVisibility','off');
end

function plot_gradient_path(ax,x,y,timeColors,isDashed,lineWidth)
    n = numel(x);
    if any(~isfinite(x)) || any(~isfinite(y))
        return;
    end
    for k = 1:(n-1)
        if isDashed && mod(floor((k-1)/4),2) == 1
            continue;
        end
        plot(ax,x(k:k+1),y(k:k+1),'-', ...
            'Color',timeColors(k,:),'LineWidth',lineWidth, ...
            'HandleVisibility','off');
    end
    scatter(ax,x(1),y(1),280,'o','MarkerFaceColor',[1 1 1], ...
        'MarkerEdgeColor',timeColors(1,:),'LineWidth',1,'HandleVisibility','off');
    scatter(ax,x(end),y(end),280,'o','MarkerFaceColor',timeColors(end,:), ...
        'MarkerEdgeColor',timeColors(end,:),'LineWidth',1,'HandleVisibility','off');
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
