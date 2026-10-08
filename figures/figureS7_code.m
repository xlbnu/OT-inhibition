%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% figure S7 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear;clc;
% Locate bundled inputs relative to this saved script, independent of pwd.
figureDir = fileparts(mfilename('fullpath'));
codeRoot = fileparts(figureDir);
dataRoot = fullfile(codeRoot,'figure_data');
addpath(fullfile(codeRoot,'utilities'));
f_file = [fullfile(codeRoot,'figure_outputs','figureS7') filesep];
if ~isfolder(f_file)
    mkdir(f_file);
end
assert(exist('ft_topoplotER','file') == 2, ...
    'Configure FieldTrip first: run setup_figure_paths(fieldtripRoot) from the package root.');

%% figure S7A-D: scalp GIFs
% (A and B) Sensor topographies and source maps of decision-related 
% conflict-specific gamma-band (A) and beta-band (B) activity. 
% (C and D) The corresponding maps of motor-related gamma-band (C) and 
% beta-band (D) activity (right minus left responses for scalp t maps).
% topoplot 
names = {'decision_gamma_band','decision_beta_band', ...
         'motor_gamma_band','motor_beta_band'};
topoplotResults = cell(size(names));
for k = 1:numel(names)
    % Motor gamma and beta use the same four top-level power variables.
    gifCfg = struct();
    gifCfg.input_dir = fullfile(dataRoot,'MEG_results','source_dynamic_gif',names{k},'topoplot');
    gifCfg.output_file = fullfile(f_file,[names{k} '_topoplot.gif']);
    topoplotResults{k} = topoplot_sliding_window_gif(names{k},gifCfg);
    assert(topoplotResults{k}.frame_count == 15,'Expected 15 scalp frames.');
end

%% figure S7A-D: source GIFs
% source map
% SPM12 is resolved from this host's MATLAB path; configure it before running.
assert(exist('spm','file') == 2, ...
    'Configure local SPM12 with setup_figure_paths(fieldtripRoot,spmRoot).');

names = {'decision_gamma_band','decision_beta_band', ...
         'motor_gamma_band','motor_beta_band'};

for k = 1:numel(names)
    result = brainnet_sliding_window_gif(names{k});
end

%% figure S7E-F
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

%% figure S7E
% source single TFR: average: gamma
int_so=all_tfr_diff1;
cfg=[];
cfg.frequency_window=[50 70];
cfg.latency_window=[-0.5 -0.3];
cfg.clim=[-3 1.5];
TFR_figure(int_so,cfg);
print(gcf,[f_file,'source_all_time-frequency_Tvalue(2-120hz)'], '-dsvg','-vector','-r600');

%% figure S7F
% decision gamma-band source reconstraction: OT-asocal
int_color = [132 121 181]./255;
time_idx=on_tfr_diff1.time>-1 & on_tfr_diff1.time<0.5;
freq_idx=on_tfr_diff1.freq>=13 & on_tfr_diff1.freq<=30;
ax1=figureX1([27 22]);
hold on
errorbound(all_tfr_diff1.time(time_idx),squeeze(mean(all_tfr_diff1.powspctrm(:,:,freq_idx,time_idx),[1 3],'omitnan')), ...
    squeeze(std(on_tfr_diff1.powspctrm(:,:,freq_idx,time_idx),[],[1 3],'omitnan'))./sqrt(size(all_tfr_diff1.powspctrm(:,:,freq_idx,time_idx),1)), ...
    'color',int_color,'linewidth',5);
plot_column_significance(all_tfr_diff1.time(time_idx),squeeze(mean(all_tfr_diff1.powspctrm(:,:,freq_idx,time_idx),3,'omitnan')), ...
    'Color',int_color,'YOffsetRatio',0.06,'LineWidth',5);
xline(0,'k--','LineWidth',1.5);yline(0,'k--','LineWidth',1.5);ylim([-0.15 0.1]);
setFixedPlotArea(ax1, [15 15],[1.8 0.8]);set(gca, 'LineWidth', 1.5,'XDir','normal'); 
hold off;xlim([-1 0.5]);
ylabel(sprintf('power(log)\n conflict-congruent'));
xlabel('time from choice made');
gcaf1(gca,35);setAxesForPPT(ax1, 44, 2.5);
print(gcf,[f_file,'decision_source_time-course_(13-30hz)_OTasocial'], '-dsvg','-vector','-r600');
%%
int_color = [132 121 181]./255;
time_idx=on_tfr_diff1.time>-1 & on_tfr_diff1.time<0.5;
freq_idx=on_tfr_diff1.freq>=60 & on_tfr_diff1.freq<=70;
tmp_tfr=on_tfr_diff1;tmp_tfr.powspctrm=(on_tfr_diff1.powspctrm)./2;
ax1=figureX1([27 22]);
hold on
errorbound(tmp_tfr.time(time_idx),squeeze(mean(tmp_tfr.powspctrm(:,:,freq_idx,time_idx),[1 3],'omitnan')), ...
    squeeze(std(tmp_tfr.powspctrm(:,:,freq_idx,time_idx),[],[1 3],'omitnan'))./sqrt(size(tmp_tfr.powspctrm(:,:,freq_idx,time_idx),1)), ...
    'color',int_color,'linewidth',5);
plot_column_significance(tmp_tfr.time(time_idx),squeeze(mean(tmp_tfr.powspctrm(:,:,freq_idx,time_idx),3,'omitnan')), ...
    'Color',int_color,'YOffsetRatio',0.04,'LineWidth',5);
xline(0,'k--','LineWidth',1.5);yline(0,'k--','LineWidth',1.5);%ylim([-0.02 0.12]);yticks([0 0.06 0.12])
setFixedPlotArea(ax1, [15 15],[1.8 0.8]);set(gca, 'LineWidth', 1.5,'XDir','normal'); 
hold off;xlim([-1 0.5]);
ylabel(sprintf('power(log)\n conflict-congruent'));
xlabel('time from choice made');
gcaf1(gca,35);setAxesForPPT(ax1, 44, 2.5);
print(gcf,[f_file,'decision_source_time-course_(60-70hz)_OTasocial'], '-dsvg','-vector','-r600');

%% figure S7G-H
load(fullfile(dataRoot,'MEG_results','motor_gamma_band_sourceTFR.mat'));
[high_all_tfr_diff,high_os_tfr_diff, ...
 high_ps_tfr_diff,high_on_tfr_diff,high_pn_tfr_diff] = ...
    getConstractPower(all_tfr_right,all_tfr_left,os_tfr_right,os_tfr_left, ...
        ps_tfr_right,ps_tfr_left,on_tfr_right,on_tfr_left,pn_tfr_right,pn_tfr_left);

cfg_ga1 = [];
cfg_ga1.keepindividual = 'yes';
all_tfr_diff2 = ft_freqgrandaverage(cfg_ga1, high_all_tfr_diff{:});
os_tfr_diff2 = ft_freqgrandaverage(cfg_ga1, high_os_tfr_diff{:});
ps_tfr_diff2 = ft_freqgrandaverage(cfg_ga1, high_ps_tfr_diff{:});
on_tfr_diff2 = ft_freqgrandaverage(cfg_ga1, high_on_tfr_diff{:});
pn_tfr_diff2 = ft_freqgrandaverage(cfg_ga1, high_pn_tfr_diff{:});
%% figure S7G
int_so=all_tfr_diff2;
cfg=[];
cfg.frequency_window=[2 100];
cfg.latency_window=[-0.8 0.5];
cfg.clim=[-5 5];
TFR_figure(int_so,cfg);
print(gcf,[f_file,'motor_source_all_time-frequency_Tvalue(2-120hz)'], '-dsvg','-vector','-r600');

%% figure S7H bottom
int_color = [132 121 181]./255;
time_idx=all_tfr_diff2.time>-1 & all_tfr_diff2.time<0.5;
freq_idx=all_tfr_diff2.freq>=13 & all_tfr_diff2.freq<=30;
ax1=figureX1([27 22]);
hold on
errorbound(all_tfr_diff2.time(time_idx),squeeze(mean(all_tfr_diff2.powspctrm(:,:,freq_idx,time_idx),[1 3],'omitnan')), ...
    squeeze(std(all_tfr_diff2.powspctrm(:,:,freq_idx,time_idx),[],[1 3],'omitnan'))./sqrt(size(all_tfr_diff2.powspctrm(:,:,freq_idx,time_idx),1)), ...
    'color',int_color,'linewidth',5);
plot_column_significance(all_tfr_diff2.time(time_idx),squeeze(mean(all_tfr_diff2.powspctrm(:,:,freq_idx,time_idx),3,'omitnan')), ...
    'Color',int_color,'YOffsetRatio',0.06,'LineWidth',5);
xline(0,'k--','LineWidth',1.5);yline(0,'k--','LineWidth',1.5);ylim([-0.14 0.20]);
setFixedPlotArea(ax1, [15 15],[1.8 0.8]);set(gca, 'LineWidth', 1.5,'XDir','normal'); 
hold off;xlim([-1 0.5]);
ylabel(sprintf('power(log)\n right-left'));
xlabel('time from choice made');
gcaf1(gca,35);setAxesForPPT(ax1, 44, 2.5);
print(gcf,[f_file,'motor_source_time-course_(13-30hz)'], '-dsvg','-vector','-r600');
% figure S7H top
int_color = [132 121 181]./255;
time_idx=all_tfr_diff2.time>-1 & all_tfr_diff2.time<0.5;
freq_idx=all_tfr_diff2.freq>=60 & all_tfr_diff2.freq<=70;
ax1=figureX1([27 22]);
hold on
errorbound(all_tfr_diff2.time(time_idx),squeeze(mean(all_tfr_diff2.powspctrm(:,:,freq_idx,time_idx),[1 3],'omitnan')), ...
    squeeze(std(all_tfr_diff2.powspctrm(:,:,freq_idx,time_idx),[],[1 3],'omitnan'))./sqrt(size(all_tfr_diff2.powspctrm(:,:,freq_idx,time_idx),1)), ...
    'color',int_color,'linewidth',5);
plot_column_significance(all_tfr_diff2.time(time_idx),squeeze(mean(all_tfr_diff2.powspctrm(:,:,freq_idx,time_idx),3,'omitnan')), ...
    'Color',int_color,'YOffsetRatio',0.2,'LineWidth',5);
xline(0,'k--','LineWidth',1.5);yline(0,'k--','LineWidth',1.5);ylim([-0.02 0.12]);yticks([0 0.06 0.12])
setFixedPlotArea(ax1, [15 15],[1.8 0.8]);set(gca, 'LineWidth', 1.5,'XDir','normal'); 
hold off;xlim([-1 0.5]);
ylabel(sprintf('power(log)\n right-left'));
xlabel('time from choice made');
gcaf1(gca,35);setAxesForPPT(ax1, 44, 2.5);
print(gcf,[f_file,'motor_source_time-course_(60-70hz)'], '-dsvg','-vector','-r600');


%% figure S7I
load(fullfile(dataRoot,'MEG_results','motor_beta_band_sourceTFR.mat'));
[~,high_os_tfr_diff, high_ps_tfr_diff,high_on_tfr_diff,high_pn_tfr_diff] = ...
    getConstractPower(all_tfr_right,all_tfr_left,os_tfr_right,os_tfr_left, ...
        ps_tfr_right,ps_tfr_left,on_tfr_right,on_tfr_left,pn_tfr_right,pn_tfr_left);

time_all=high_os_tfr_diff{1}.time;
freq_all=high_os_tfr_diff{1}.freq;
time_index=(round(time_all,3)>=-0.2) & (round(time_all,3)<=0);
freq_index=(round(freq_all)>=13) & (round(freq_all)<=30);
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
box off;ylabel(sprintf('power(log)\n right-left'));
set(gca,'XTick',x_pos,'XTickLabel',{'OT','PL','OT','PL'}); gcaf1(gca,35);
setFixedPlotArea(ax1, [15 15]);setAxesForPPT(ax1, 35, 2.5);
cutBackColor(gca,gcf);print(gcf,[f_file,'motor_beta_tfr_interaction(13-30hz)'], '-dsvg','-vector','-r600');

%% figure S7J
load(fullfile(dataRoot,'MEG_results','motor_gamma_band_sourceTFR.mat'));
[high_all_tfr_diff,high_os_tfr_diff, ...
 high_ps_tfr_diff,high_on_tfr_diff,high_pn_tfr_diff] = ...
    getConstractPower(all_tfr_right,all_tfr_left,os_tfr_right,os_tfr_left, ...
        ps_tfr_right,ps_tfr_left,on_tfr_right,on_tfr_left,pn_tfr_right,pn_tfr_left);

time_all=high_os_tfr_diff{1}.time;
freq_all=high_os_tfr_diff{1}.freq;
time_index=(round(time_all,3)>=-0.1) & (round(time_all,3)<=0.1);
freq_index=(round(freq_all)>=40) & (round(freq_all)<=70);
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
box off;ylabel(sprintf('power(log)\n right-left'));
set(gca,'XTick',x_pos,'XTickLabel',{'OT','PL','OT','PL'}); gcaf1(gca,35);
setFixedPlotArea(ax1, [15 15]);setAxesForPPT(ax1, 35, 2.5);
cutBackColor(gca,gcf);print(gcf,[f_file,'motor_gamma_tfr_interaction(40-70hz)'], '-dsvg','-vector','-r600');



%% figure S7K
% motor & decision latency
choice_gamma=all_tfr_diff1;motor_gamma=all_tfr_diff2;
time_index1=choice_gamma.time>=-1 & choice_gamma.time<=0;
time_index2=motor_gamma.time>=-1 & motor_gamma.time<=0.5;

x_1=choice_gamma.time(time_index1);x_2=choice_gamma.time(time_index1);
x_3=motor_gamma.time(time_index2);x_4=motor_gamma.time(time_index2);

freq_index2=round(choice_gamma.freq)>=60 & round(choice_gamma.freq)<=70;
freq_index3=round(motor_gamma.freq)>=13 & round(motor_gamma.freq)<=30;
freq_index1=round(choice_gamma.freq)>=13 & round(choice_gamma.freq)<=30;
freq_index4=round(motor_gamma.freq)>=60 & round(motor_gamma.freq)<=70;

mgamma_peak=nan(19,1);mgamma_latency=nan(19,1);
mbeta_peak=nan(19,1);mbeta_latency=nan(19,1);
dgamma_peak=nan(19,1);dgamma_latency=nan(19,1);
dbeta_peak=nan(19,1);dbeta_latency=nan(19,1);
for i=1:19
data1=squeeze(mean(choice_gamma.powspctrm(:,:,freq_index2,time_index1),3,'omitnan'));
[dgamma_peak(i,1),tmp1]=max((data1(i,:)));
dgamma_latency(i,1)=x_2(tmp1);

data2=squeeze(mean(choice_gamma.powspctrm(:,:,freq_index1,time_index1),3,'omitnan'));
[dbeta_peak(i,1),tmp1]=min(data2(i,:));
dbeta_latency(i,1)=x_1(tmp1);

data3=squeeze(mean(motor_gamma.powspctrm(:,:,freq_index4,time_index2),3,'omitnan'));
[mgamma_peak(i,1),tmp1]=max(data3(i,:));
mgamma_latency(i,1)=x_4(tmp1);


data4=squeeze(mean(motor_gamma.powspctrm(:,:,freq_index3,time_index2),3,'omitnan'));
[mbeta_peak(i,1),tmp1]=min(data4(i,:));
mbeta_latency(i,1)=x_3(tmp1);
end

%
ax1=figureX1([60 20]);
xx1=[dbeta_latency,dgamma_latency,mbeta_latency,mgamma_latency];
x_pos=figureTmp1_nonGroup(xx1);
box off;ylabel(sprintf('time (peak)'));
ax = gca;
set(ax, 'XTick', x_pos, 'XTickLabel', {'decision-beta','decision-gamma','motor-beta','motor-gamma'},'XTickLabelRotation',0);
gcaf1(gca,35);
setFixedPlotArea(ax1, [52.5 15]);setAxesForPPT(ax1, 40, 2.5);
print(gcf,[f_file,'choice&motor source lantency'], '-dsvg','-vector','-r600');
%































%%
%%%%%%%%%%%%%%%%%%%%%%%%%%% function %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
