%% Paired source F-test on cerebral cortex and interpolation
%
% Before running this script, provide two paired source cell arrays:
%
%   source_condition1 = ...;   % N subjects x 1 cell
%   source_condition2 = ...;   % N subjects x 1 cell
%
% The current WTA source program stores the source value in .effect, so the
% parameter may be 'effect'. Select the actual parameter for each section:
% the SigmaV and DICS branches below load external inputs and export results.
% Configure their load/save paths before running them.

%% coefficient(ΣV) source cluster-permutation F test
load(['F:\xianliang\exp_data\OT_N2_MEG\MEG_source_408IndependDay\beta_plus\' ...
    'SigmaV_betacov_F(T-0.50_-0.30)_T(-0.50_-0.30)_FI(0.00_0.20)_TI(0.00_0.20)_u1\SigmaV_vectorMagnitude_and_SumSq_source.mat'])
for i=1:19
subject_result.by_measure.sumsq.target_unitnoisegain{i}.pow=log(subject_result.by_measure.sumsq.target_unitnoisegain{i}.pow);
subject_result.by_measure.sumsq.iti_unitnoisegain{i}.pow=log(subject_result.by_measure.sumsq.iti_unitnoisegain{i}.pow);
end
source_condition1 = subject_result.by_measure.sumsq.target_unitnoisegain;
source_condition2 =subject_result.by_measure.sumsq.iti_unitnoisegain;

cfg_source_stat = struct();
cfg_source_stat.parameter = 'pow';
cfg_source_stat.correctm = 'cluster';
cfg_source_stat.clusteralpha = 0.05;
cfg_source_stat.alpha = 0.05;
cfg_source_stat.numrandomization = 5000;
[stat_res,template_grid,final_inside,atlas,atlas_grid]=clusterFtest(cfg_source_stat,source_condition1,source_condition2);
%
idx_cluster_all=getTargetCluster(1,stat_res,atlas,atlas_grid);

% Export the NIfTI statistic for the first significant cluster after cluster correction
ExportNIfTI(['F:\xianliang\exp_data\OT_N2_MEG\MEG_source_408IndependDay\beta_plus\' ...
    'SigmaV_betacov_F(T-0.50_-0.30)_T(-0.50_-0.30)_FI(0.00_0.20)_TI(0.00_0.20)_u1\all_gamma_source_cluster']);
% Export the NIfTI statistic for oxytocin receptor expression map conjunction cluster after cluster correction
out_put=['F:\xianliang\exp_data\OT_N2_MEG\MEG_source_408IndependDay\beta_plus\' ...
    'SigmaV_betacov_F(T-0.50_-0.30)_T(-0.50_-0.30)_FI(0.00_0.20)_TI(0.00_0.20)_u1'];
export_significant_cluster1_oxtr_to_nifti(stat_res,fullfile(out_put,'conjunction_'));
%% gamma-band source cluster-permutation F test
load('F:\xianliang\exp_data\OT_N2_MEG\MEG_source_408IndependDay\unitGain\Source_408DICS_rank3_noDayUnitGain_unitWhiteNoiseGain_30-100Hz_(-0.5_-0.3)\all_source.mat')
for i=1:19
all_source_conf{i}.avg.pow=log(all_source_conf{i}.avg.pow);
all_source_cong{i}.avg.pow=log(all_source_cong{i}.avg.pow);
end
source_condition1 = all_source_conf;
source_condition2 =all_source_cong;

cfg_source_stat = struct();
cfg_source_stat.parameter = 'avg.pow';
cfg_source_stat.correctm = 'cluster';
cfg_source_stat.clusteralpha = 0.01;
cfg_source_stat.alpha = 0.05;
cfg_source_stat.numrandomization = 5000;
[stat_res,template_grid,final_inside,atlas,atlas_grid]=clusterFtest(cfg_source_stat,source_condition1,source_condition2);
%
idx_cluster_all=getTargetCluster(1,stat_res,atlas,atlas_grid);

% Export the NIfTI statistic for the first significant cluster after cluster correction
ExportNIfTI(['F:\xianliang\exp_data\OT_N2_MEG\MEG_source_408IndependDay\unitGain\' ...
    'Source_408DICS_rank3_noDayUnitGain_unitWhiteNoiseGain_30-100Hz_(-0.5_-0.3)\all_gamma_source_cluster']);

% Export the NIfTI statistic for oxytocin receptor expression map conjunction cluster after cluster correction
out_put=['F:\xianliang\exp_data\OT_N2_MEG\MEG_source_408IndependDay\unitGain\' ...
    'Source_408DICS_rank3_noDayUnitGain_unitWhiteNoiseGain_30-100Hz_(-0.5_-0.3)'];
export_significant_cluster1_oxtr_to_nifti(stat_res,fullfile(out_put,'conjunction_'));

%%  Interpolate the F statistic and corrected binary mask
ft_root = fileparts(which('ft_defaults'));
template_mri = ft_read_mri(fullfile(ft_root,'template','anatomy', ...
    'single_subj_T1_1mm.nii'));
template_mri = ft_convert_units(template_mri,'mm');

stat_clean = struct();
stat_clean.pos = template_grid.pos;
stat_clean.dim = template_grid.dim;
stat_clean.inside = final_inside;
stat_clean.unit = template_grid.unit;
stat_clean.stat = double(stat_res.stat(:));
stat_clean.stat(~final_inside | ~isfinite(stat_clean.stat)) = 0;

cfg_interp = [];
cfg_interp.parameter = 'stat';
cfg_interp.interpmethod = 'linear';
stat_interp = ft_sourceinterpolate(cfg_interp,stat_clean,template_mri);

mask_clean = rmfield(stat_clean,'stat');
mask_clean.mask = double(stat_res.mask);
cfg_interp = [];
cfg_interp.parameter = 'mask';
cfg_interp.interpmethod = 'linear';
mask_interp = ft_sourceinterpolate(cfg_interp,mask_clean,template_mri);
stat_interp.mask = isfinite(mask_interp.mask) & mask_interp.mask>=0.5;

disp('F-test and interpolation completed: stat_res, stat_interp.');


%





















%%
function idx_cluster_all=getTargetCluster(target_cluster_id,stat_res,atlas,atlas_grid)
disp('▶ 正在生成 Cluster 空间分布详细清单...');

% 1. Specify the target to query Cluster identifier
% target_cluster_id = 1; 

% 2. Extract this Cluster all included grid-point indices
idx_cluster_all = find(stat_res.posclusterslabelmat == target_cluster_id);

if isempty(idx_cluster_all)
    error('未找到指定的 Cluster ID！');
end

% 3. Get these grid points' AAL corresponding atlas tissue IDs
cluster_tissue_indices = atlas_grid.tissue(idx_cluster_all);

% 4. Exclude points not covered by the atlas"Unknown regions"(tissue ID is 0 points, Typically in ventricles or skull boundaries)
valid_tissue_indices = cluster_tissue_indices(cluster_tissue_indices > 0);
num_unlabeled = sum(cluster_tissue_indices == 0);

% 5. Find this Cluster all involved distinct brain-region IDs
unique_tissue_ids = unique(valid_tissue_indices);

% 6. Count grid points in each brain region
region_counts = zeros(length(unique_tissue_ids), 1);
region_names = cell(length(unique_tissue_ids), 1);

for i = 1:length(unique_tissue_ids)
    tissue_id = unique_tissue_ids(i);
    region_counts(i) = sum(valid_tissue_indices == tissue_id);
    region_names{i} = atlas.tissuelabel{tissue_id};
end

% 7. 🌟 Core optimization: Sort by descending number of included grid points
[sorted_counts, sort_idx] = sort(region_counts, 'descend');
sorted_names = region_names(sort_idx);

% Compute the total valid grid count, For percentage calculation
total_valid_voxels = sum(sorted_counts);

% 8. Print a detailed report
fprintf('\n===================================================\n');
fprintf('🎯 【Cluster %d 空间分布解剖报告】\n', target_cluster_id);
fprintf('   总网格点数 : %d (其中 %d 个点不在 AAL 标签范围内)\n', length(idx_cluster_all), num_unlabeled);
fprintf('---------------------------------------------------\n');
fprintf('%-25s | %-10s | %-10s\n', '解剖脑区名称 (AAL)', '网格数量', '所占比例');
fprintf('---------------------------------------------------\n');

for i = 1:length(sorted_names)
    percentage = (sorted_counts(i) / total_valid_voxels) * 100;
    
    % Print only regions exceeding the proportion 1% or containing at least 2 grid points, Filter very small peripheral extensions
    if sorted_counts(i) >= 2 || percentage >= 1.0
        fprintf('%-27s | %6d     | %6.1f%%\n', sorted_names{i}, sorted_counts(i), percentage);
    end
end
fprintf('===================================================\n\n');
end





function [stat_res,template_grid,final_inside,atlas,atlas_grid]=clusterFtest(cfg_source_stat,source_condition1,source_condition2)
%% 1. Statistical parameters
if ~exist('cfg_source_stat','var') || isempty(cfg_source_stat)
    cfg_source_stat = struct();
end

if ~isfield(cfg_source_stat,'parameter')
    cfg_source_stat.parameter = 'effect';  % effect | avg.pow
end
if ~isfield(cfg_source_stat,'correctm')
    cfg_source_stat.correctm = 'cluster';  % cluster | no
end
if ~isfield(cfg_source_stat,'clusteralpha')
    cfg_source_stat.clusteralpha = 0.05;
end
if ~isfield(cfg_source_stat,'alpha')
    cfg_source_stat.alpha = 0.05;
end
if ~isfield(cfg_source_stat,'numrandomization')
    cfg_source_stat.numrandomization = 5000;
end
if ~isfield(cfg_source_stat,'minnbchan')
    cfg_source_stat.minnbchan = 0;
end

assert(exist('source_condition1','var')==1, ...
    'Please provide source_condition1 in the workspace.');
assert(exist('source_condition2','var')==1, ...
    'Please provide source_condition2 in the workspace.');
assert(iscell(source_condition1) && iscell(source_condition2), ...
    'Both source inputs must be cell arrays.');
assert(numel(source_condition1)==numel(source_condition2), ...
    'The two conditions contain different subject counts.');
assert(all(~cellfun(@isempty,source_condition1)) && ...
       all(~cellfun(@isempty,source_condition2)), ...
    'One or more paired source entries are empty.');

num_subjects = numel(source_condition1);
assert(num_subjects>=2,'At least two paired subjects are required.');
assert(exist('ft_defaults','file')==2, ...
    'FieldTrip must already be on the MATLAB path.');
ft_defaults;

% 2. Generate the cerebral-cortex 5-mm grid
% Locate FieldTrip templates from the current MATLAB path. No hard-coded
% FieldTrip directory is required.
ft_root = fileparts(which('ft_defaults'));
atlas = ft_read_atlas(fullfile(ft_root,'template','atlas','aal', ...
    'ROI_MNI_V4.nii'));

grid_file = load(fullfile(ft_root,'template','sourcemodel', ...
    'standard_sourcemodel3d5mm.mat'));
grid_names = fieldnames(grid_file);
template_grid = grid_file.(grid_names{1});
template_grid = ft_convert_units(template_grid,'mm');
n_grid = size(template_grid.pos,1);

cfg_atlas = [];
cfg_atlas.interpmethod = 'nearest';
cfg_atlas.parameter = 'tissue';
atlas_grid = ft_sourceinterpolate(cfg_atlas,atlas,template_grid);

% Exclude cerebellum, brainstem and deep nuclei. Lingual gyrus is cortex
% and is therefore not excluded.
exclude_keywords = { ...
    'Cerebelum','Cerebellum','Vermis','Brainstem','Brain_Stem', ...
    'Thalamus','Caudate','Putamen','Pallidum','Amygdala', ...
    'Hippocampus','Accumbens','VentralDC','Substantia_Nigra', ...
    'Red_Nucleus','Olfactory','Lingual'};%

cortical_mask = false(n_grid,1);
for i_roi = 1:numel(atlas.tissuelabel)
    roi_name = atlas.tissuelabel{i_roi};
    exclude_roi = false;
    for i_key = 1:numel(exclude_keywords)
        if contains(roi_name,exclude_keywords{i_key},'IgnoreCase',true)
            exclude_roi = true;
            break;
        end
    end
    if ~exclude_roi
        cortical_mask(double(atlas_grid.tissue(:))==i_roi) = true;
    end
end

if islogical(template_grid.inside)
    template_inside = template_grid.inside(:);
else
    template_inside = false(n_grid,1);
    template_inside(double(template_grid.inside(:))) = true;
end

% 3. Extract values and retain only common valid cortical grid points
condition1_value = nan(num_subjects,n_grid);
condition2_value = nan(num_subjects,n_grid);
condition1_inside = false(num_subjects,n_grid);
condition2_inside = false(num_subjects,n_grid);

for i_sub = 1:num_subjects
    condition1_value(i_sub,:) = get_source_parameter( ...
        source_condition1{i_sub},cfg_source_stat.parameter,n_grid);
    condition2_value(i_sub,:) = get_source_parameter( ...
        source_condition2{i_sub},cfg_source_stat.parameter,n_grid);
    condition1_inside(i_sub,:) = get_source_inside( ...
        source_condition1{i_sub},n_grid);
    condition2_inside(i_sub,:) = get_source_inside( ...
        source_condition2{i_sub},n_grid);
end

final_inside = template_inside & cortical_mask & ...
    all(condition1_inside,1)' & all(condition2_inside,1)' & ...
    all(isfinite(condition1_value),1)' & ...
    all(isfinite(condition2_value),1)';
assert(any(final_inside),'No common valid cortical grid points remain.');

cortical_grid = template_grid;
cortical_grid.inside = final_inside;
fprintf('Subjects: %d; tested cortical grid points: %d\n', ...
    num_subjects,sum(final_inside));

% Repackage both conditions into one uniform FieldTrip source format.
condition1_ft = cell(num_subjects,1);
condition2_ft = cell(num_subjects,1);
for i_sub = 1:num_subjects
    condition1_ft{i_sub} = template_grid;
    condition1_ft{i_sub}.inside = final_inside;
    condition1_ft{i_sub}.pow = condition1_value(i_sub,:)';
    condition1_ft{i_sub}.pow(~final_inside) = NaN;
    condition1_ft{i_sub}.dimord = 'pos';

    condition2_ft{i_sub} = template_grid;
    condition2_ft{i_sub}.inside = final_inside;
    condition2_ft{i_sub}.pow = condition2_value(i_sub,:)';
    condition2_ft{i_sub}.pow(~final_inside) = NaN;
    condition2_ft{i_sub}.dimord = 'pos';
end

% 4. Direct paired F-test between the two source conditions
% For exactly two paired conditions, the following are equivalent when the
% same subjects, preprocessing, mask and permutation scheme are used:
%
%   F(condition1 versus condition2)
%   F((condition1-condition2) versus zero)
%   paired_t(condition1 versus condition2)^2
%
% Direct comparison is used here. F is nonnegative and does not indicate
% which condition is larger; inspect the condition1-condition2 mean if
% direction is needed.
cfg_stat = [];
cfg_stat.dim = template_grid.dim;
cfg_stat.method = 'montecarlo';
cfg_stat.statistic = 'ft_statfun_depsamplesFunivariate';
cfg_stat.parameter = 'pow';
cfg_stat.correctm = cfg_source_stat.correctm;
cfg_stat.tail = 1;
cfg_stat.alpha = cfg_source_stat.alpha;
cfg_stat.numrandomization = cfg_source_stat.numrandomization;

if strcmpi(cfg_source_stat.correctm,'cluster')
    cfg_stat.clusteralpha = cfg_source_stat.clusteralpha;
    cfg_stat.clusterstatistic = 'maxsum';
    cfg_stat.clustertail = 1;
    cfg_stat.minnbchan = cfg_source_stat.minnbchan;
end

cfg_stat.design = zeros(2,2*num_subjects);
cfg_stat.design(1,:) = [1:num_subjects,1:num_subjects];
cfg_stat.design(2,:) = [ones(1,num_subjects),2*ones(1,num_subjects)];
cfg_stat.uvar = 1;
cfg_stat.ivar = 2;

stat_res = ft_sourcestatistics(cfg_stat, ...
    condition1_ft{:},condition2_ft{:});

% Sanitize the corrected mask; avoids NaN-to-logical conversion errors.
if isfield(stat_res,'mask')
    stat_res.mask = isfinite(double(stat_res.mask(:))) & ...
        double(stat_res.mask(:))~=0 & final_inside;
else
    stat_res.mask = false(n_grid,1);
end
end



%% Local helpers in the same script file
function value = get_source_parameter(source,parameter,n_grid)
parts = strsplit(char(string(parameter)),'.');
value = source;
for i = 1:numel(parts)
    assert(isstruct(value) && isfield(value,parts{i}), ...
        'Source structure lacks parameter %s.',parameter);
    value = value.(parts{i});
end
assert(isnumeric(value) && numel(value)==n_grid, ...
    'Source parameter %s must contain %d grid values.',parameter,n_grid);
value = double(value(:))';
end


function inside = get_source_inside(source,n_grid)
if ~isfield(source,'inside') || isempty(source.inside)
    inside = true(1,n_grid);
    return;
end
x = source.inside;
if islogical(x)
    assert(numel(x)==n_grid,'Source inside mask has the wrong length.');
    inside = x(:)';
elseif numel(x)==n_grid && ...
        all(ismember(unique(double(x(isfinite(x)))),[0 1]))
    inside = (isfinite(double(x(:))) & double(x(:))~=0)';
else
    inside = false(1,n_grid);
    idx = double(x(:));
    idx = idx(isfinite(idx) & idx==round(idx) & idx>=1 & idx<=n_grid);
    inside(idx) = true;
end
end

























%%
% RUN_WTA_SOURCE_F_CLUSTER_STATISTICS
%
% Group-level cortical-mask cluster permutation F test for the current
% time-domain WTA source-localization output.
%
% Default tested quantity:
%   [(beta_winner-beta_loser)_target -
%    (beta_winner-beta_loser)_ITI] versus zero.
%
% F is non-negative and does not encode direction. Inspect
% stat_res.mean_effect (and cluster_summary.mean_effect) to determine
% whether a significant cluster is positive or negative.

cfg_main = struct();
cfg_main.ft_path = 'D:\software\matlab_toolbox\fieldtrip-20251218';
cfg_main.group_file = ['F:\xianliang\exp_data\OT_N2_MEG\MEG_source_408IndependDay\WTA_beta_rawcov_jointLCMV_T(-0.5_-0.3)_I(0.2)delay0.15\' ...
    'GROUP_WTA_beta_rawcov_jointLCMV.mat'];
cfg_main.template_grid_file = fullfile(cfg_main.ft_path,'template', ...
    'sourcemodel','standard_sourcemodel3d5mm.mat');
cfg_main.aal_file = fullfile(cfg_main.ft_path,'template','atlas', ...
    'aal','ROI_MNI_V4.nii');
cfg_main.template_mri_file = fullfile(cfg_main.ft_path,'template', ...
    'anatomy','single_subj_T1_1mm.nii');
cfg_main.metric = 'signed_delta';
cfg_main.clusteralpha = 0.01;
cfg_main.alpha = 0.05;
cfg_main.numrandomization = 1000;
cfg_main.minnbchan = 0;
cfg_main.randomseed = 20260717;
cfg_main.interpolate = true;

cfg_main.metric = char(lower(string(cfg_main.metric)));
valid_metric = {'signed_delta','magnitude_delta', ...
    'target_signed','target_magnitude'};
assert(ismember(cfg_main.metric,valid_metric), ...
    'cfg.metric must be one of: %s.',strjoin(valid_metric,', '));


addpath(cfg_main.ft_path);
ft_defaults;
rng(cfg_main.randomseed,'twister');


%% Load current source-localization results
S = load(cfg_main.group_file,'all_source_target','all_source_iti', ...
    'all_source_delta','subject_ids');
required = {'all_source_target','all_source_iti','all_source_delta'};
for i = 1:numel(required)
    assert(isfield(S,required{i}),'Group file lacks %s.',required{i});
end

n_sub = numel(S.all_source_delta);
assert(n_sub >= 2,'At least two subjects are required.');
assert(numel(S.all_source_target)==n_sub && ...
    numel(S.all_source_iti)==n_sub,'Source cell counts differ.');
if isfield(S,'subject_ids') && numel(S.subject_ids)==n_sub
    subject_ids = S.subject_ids;
else
    subject_ids = arrayfun(@(x)sprintf('sub%02d',x), ...
        (1:n_sub)','UniformOutput',false);
end

%% Load the exact 5-mm template grid and standardize geometry
G = load(cfg_main.template_grid_file);
assert(isfield(G,'sourcemodel'),'Template file lacks sourcemodel.');
template_grid = ft_convert_units(G.sourcemodel,'mm');
n_grid = size(template_grid.pos,1);
template_inside = inside_to_logical(template_grid.inside,n_grid);

%% Construct a pure cerebral-cortex mask from AAL
atlas = ft_read_atlas(cfg_main.aal_file);
atlas = ft_convert_units(atlas,'mm');

cfg_interp = [];
cfg_interp.interpmethod = 'nearest';
cfg_interp.parameter = 'tissue';
atlas_grid = ft_sourceinterpolate(cfg_interp,atlas,template_grid);

% AAL116 has no explicit brainstem label. Unlabelled template points are
% never added to this mask, so brainstem/white matter are excluded.
exclude_keywords = lower([ ...
    "Cerebellum","Vermis", ...
    "Thalamus","Caudate","Putamen","Pallidum", ...
    "Amygdala","Hippocampus"]);
% exclude_keywords = lower(["Cerebellum"]);
cortical_mask = false(n_grid,1);
kept_labels = {};
excluded_labels = {};
for i_roi = 1:numel(atlas.tissuelabel)
    roi_name = string(atlas.tissuelabel{i_roi});
    is_excluded = any(contains(lower(roi_name),exclude_keywords));
    if is_excluded
        excluded_labels{end+1,1} = char(roi_name); 
    else
        cortical_mask(atlas_grid.tissue(:)==i_roi) = true;
        kept_labels{end+1,1} = char(roi_name); 
    end
end
cortical_mask = cortical_mask & template_inside;

fprintf('\n============================================================\n');
fprintf('Cortical-mask WTA source F test\n');
fprintf('Metric              : %s\n',cfg_main.metric);
fprintf('Template inside     : %d points\n',sum(template_inside));
fprintf('Pure cortical mask  : %d points\n',sum(cortical_mask));
fprintf('Excluded labels     : %d\n',numel(excluded_labels));
fprintf('============================================================\n');

%% Extract the subject-level statistic and construct common valid mask
effect_matrix = nan(n_grid,n_sub);
source_effect_all = cell(n_sub,1);
source_zero_all = cell(n_sub,1);
common_inside = cortical_mask;

for i_sub = 1:n_sub
    target = standardize_source(S.all_source_target{i_sub}, ...
        template_grid,sprintf('target subject %d',i_sub));
    iti = standardize_source(S.all_source_iti{i_sub}, ...
        template_grid,sprintf('ITI subject %d',i_sub));
    delta = standardize_source(S.all_source_delta{i_sub}, ...
        template_grid,sprintf('delta subject %d',i_sub));

    effect = select_effect(target,iti,delta,cfg_main.metric);
    effect = double(effect(:));
    assert(numel(effect)==n_grid);

    valid_source = inside_to_logical(target.inside,n_grid) & ...
        inside_to_logical(iti.inside,n_grid) & ...
        inside_to_logical(delta.inside,n_grid) & isfinite(effect);
    common_inside = common_inside & valid_source;
    effect_matrix(:,i_sub) = effect;
end

assert(any(common_inside),'No common cortical source points remain.');
fprintf('Common valid cortex : %d points\n',sum(common_inside));

% Every subject must use exactly the same statistical search space.
for i_sub = 1:n_sub
    effect = effect_matrix(:,i_sub);
    effect(~common_inside) = NaN;
    source_effect_all{i_sub} = make_ft_source( ...
        template_grid,common_inside,effect);
    source_zero_all{i_sub} = make_ft_source( ...
        template_grid,common_inside,zeros(n_grid,1));
end

%% Paired two-level F test: subject effect map versus zero map
cfg_stat = [];
cfg_stat.method = 'montecarlo';
cfg_stat.statistic = 'ft_statfun_depsamplesFunivariate';
cfg_stat.parameter = 'avg.pow';
cfg_stat.correctm = 'cluster';
cfg_stat.clusteralpha = cfg_main.clusteralpha;
cfg_stat.clusterstatistic = 'maxsum';
cfg_stat.minnbchan = cfg_main.minnbchan;
cfg_stat.tail = 1;          % F >= 0
cfg_stat.clustertail = 1;   % F >= 0
cfg_stat.alpha = cfg_main.alpha;
cfg_stat.numrandomization = cfg_main.numrandomization;
cfg_stat.dim = template_grid.dim;

cfg_stat.design = zeros(2,2*n_sub);
cfg_stat.design(1,:) = [1:n_sub 1:n_sub];
cfg_stat.design(2,:) = [ones(1,n_sub) 2*ones(1,n_sub)];
cfg_stat.uvar = 1;
cfg_stat.ivar = 2;

stat_res = ft_sourcestatistics(cfg_stat, ...
    source_effect_all{:},source_zero_all{:});

% F values do not carry sign. Preserve the group mean effect and Cohen dz.
mean_effect = mean(effect_matrix,2,'omitnan');
sd_effect = std(effect_matrix,0,2,'omitnan');
cohen_dz = mean_effect./sd_effect;
cohen_dz(sd_effect<=eps | ~isfinite(cohen_dz)) = NaN;
mean_effect(~common_inside) = NaN;
cohen_dz(~common_inside) = NaN;
stat_res.mean_effect = mean_effect;
stat_res.cohen_dz = cohen_dz;
stat_res.mean_effectdimord = 'pos';
stat_res.cohen_dzdimord = 'pos';

cluster_summary = summarize_F_clusters(stat_res,mean_effect, ...
    template_grid.pos,cfg_main.alpha);
disp(cluster_summary);

%% Interpolate F, direction and corrected mask for visualization
stat_interp = [];
if cfg_main.interpolate
    template_mri = ft_read_mri(cfg_main.template_mri_file);
    template_mri = ft_convert_units(template_mri,'mm');

    stat_clean = struct();
    stat_clean.pos = stat_res.pos;
    stat_clean.dim = stat_res.dim;
    stat_clean.inside = stat_res.inside;
    stat_clean.unit = stat_res.unit;
    stat_clean.stat = stat_res.stat;
    stat_clean.mean_effect = stat_res.mean_effect;

    cfg_i = [];
    cfg_i.parameter = {'stat','mean_effect'};
    cfg_i.interpmethod = 'linear';
    stat_interp = ft_sourceinterpolate(cfg_i,stat_clean,template_mri);

    if isfield(stat_res,'mask')
        corrected_mask = isfinite(stat_res.mask) & stat_res.mask ~= 0;
    else
        corrected_mask = false(size(stat_res.stat));
    end
    mask_clean = stat_clean;
    mask_clean = rmfield(mask_clean,{'stat','mean_effect'});
    mask_clean.mask = double(corrected_mask);
    cfg_i = [];
    cfg_i.parameter = 'mask';
    cfg_i.interpmethod = 'nearest';
    mask_interp = ft_sourceinterpolate(cfg_i,mask_clean,template_mri);
    stat_interp.mask = logical(mask_interp.mask>0.5);
end

%%

%%
custom_cmap = zeros(256, 3);
custom_cmap(:, 1) = 1;                              % R Always 1 (Fully red)
custom_cmap(:, 2) = linspace(1, 0, 256)';           % G From 1 Decrease to 0 (Red+Green=Yellow, Decrease to 0 Become pure red)
custom_cmap(:, 3) = 0;   
custom_cmap=custom_cmap(end:-1:1,:);
%%
figure('Name', 'Group Level: True Cortical Activation (Ortho)', 'Color', 'w');

cfg_plot                = [];
cfg_plot.method         = 'ortho';          
cfg_plot.funparameter   = 'stat';           
cfg_plot.maskparameter  = 'mask';           
cfg_plot.funcolormap    = custom_cmap;
cfg_plot.funcolorlim    = 'maxabs';         
stat_interp.coordsys = 'mni';
% ==========================================
% 💡 [Core fix]: Explicitly set the opacity range, Avoid ALim Crash
cfg_plot.opacitylim     = [0 1];       
cfg_plot.funcolorlim = [0 20];

% cfg_plot.opacitymap   = 'rampup'; % ⚠️ Comment out this line, Because mask gradient opacity is no longer needed
% ==========================================
cfg_plot.atlas          = atlas;       
% Strongly recommended to add this: Automatically move the crosshair to T the maximum value(most significant)brain region!
% cfg_plot.location       = 'max';            
% cfg_plot.roi = {'Cingulum Ant L','Frontal Sup Medial L'};
ft_sourceplot(cfg_plot, stat_interp);
%% === Combined with AAL interactive atlas visualization (Version fixing missing-state errors) ===
figure('Name', 'Interaction Effect with Atlas', 'Color', 'w');

stat_interp.coordsys = 'mni';

% 2. Configure plotting parameters
cfg_plot                = [];
cfg_plot.method         = 'ortho';          
cfg_plot.funparameter   = 'stat';           
cfg_plot.maskparameter  = 'inside';           % Use gradient mapping to visualize trends
cfg_plot.opacitymap     = 'vdown';          
cfg_plot.opacitylim     = 'maxabs';         
cfg_plot.funcolormap    = 'jet';
cfg_plot.funcolorlim    = [0 10];           
% cfg_plot.location       = 'min';            % Automatically select the strongest signal point
cfg_plot.location      = 'center'; % Automatically locate the energy center
% Add the atlas to the renderer
cfg_plot.atlas          = atlas;            

ft_sourceplot(cfg_plot, stat_interp);

%% === Plot a hard-thresholded (|T| > 2.1) 's Conflict effect map (Ortho view) ===

% 💡 If the data structure name is not stat_interp, use the actual name
stat_unthresholded = stat_interp; 

% 💡 Set the desired fixed hard threshold here (For example, uncorrected p=0.05 corresponding T value)
t_critical_threshold = 2; 

fprintf('\n▶ 正在生成硬阈值 (|T| > %.1f) 的二值掩码...\n', t_critical_threshold);

% ==========================================================
% 💡 [Core step one]: in MATLAB Generate a binary mask in memory
% ==========================================================
% Use MATLAB logical operator.
% Create a field named 'hard_mask' in the data structure.
% |T| > 2.1 locations are assigned 1 (true, opaque); Assign the remaining whole brain 0 (false, fully transparent).
stat_unthresholded.hard_mask = abs(stat_unthresholded.stat) > t_critical_threshold;


% ==========================================================
% 💡 [Core step two]: configuration ft_sourceplot Use this binary mask
% ==========================================================
figure('Name', sprintf('Hard Thresholded |T| > %.1f (Uncorrected)', t_critical_threshold), 'Color', 'w');

cfg_plot                = [];
cfg_plot.method         = 'ortho';          
cfg_plot.funparameter   = 'stat';           % Colors still represent the actual T value gradient, For example, red-Yellow-Green-Blue

% 🌟 [Final opacity settings]: Implement"on/off"hard-threshold display
cfg_plot.maskparameter  = 'hard_mask';      % Explicitly use the newly created hard binary mask
cfg_plot.opacitymap     = 'rampup';         % 0 Map to full transparency, 1 Map to full opacity
cfg_plot.opacitylim     = [0 1];            % Fix the binary mapping range

cfg_plot.funcolormap    = 'autumn';
% cfg_plot.funcolorlim    = [0 15];%[0 10];           % Recommended fixed color-bar range, Make the gradient more pronounced

% Add automatic atlas localization labels(Avoid coordsys error-fixing version)
stat_unthresholded.coordsys = 'mni';
% ft_path = 'D:\software\matlab_toolbox\fieldtrip-20251218';
% cfg_plot.atlas = ft_read_atlas(fullfile(ft_path, 'template', 'atlas', 'aal', 'ROI_MNI_V4.nii'));
cfg_plot.atlas          = atlas;   
% cfg_plot.location = 'max'; % Automatically select the strongest signal point
cfg.location      = 'center'; % Automatically locate the energy center
% Plotting
ft_sourceplot(cfg_plot, stat_unthresholded);

disp('  ✔ 硬阈值图像已生成！');


























function cfg = set_defaults(cfg)
defaults = struct();
defaults.ft_path = 'D:\software\matlab_toolbox\fieldtrip-20251218';
defaults.group_file = ['F:\xianliang\exp_data\OT_N2_MEG\' ...
    'MEG_source_408IndependDay\WTA_beta_rawcov_jointLCMV\' ...
    'GROUP_WTA_beta_rawcov_jointLCMV.mat'];
defaults.output_dir = ['F:\xianliang\exp_data\OT_N2_MEG\' ...
    'MEG_source_408IndependDay\WTA_beta_rawcov_jointLCMV\statistics'];
defaults.template_grid_file = fullfile(defaults.ft_path,'template', ...
    'sourcemodel','standard_sourcemodel3d5mm.mat');
defaults.aal_file = fullfile(defaults.ft_path,'template','atlas', ...
    'aal','ROI_MNI_V4.nii');
defaults.template_mri_file = fullfile(defaults.ft_path,'template', ...
    'anatomy','single_subj_T1_1mm.nii');
defaults.metric = 'signed_delta';
defaults.clusteralpha = 0.01;
defaults.alpha = 0.05;
defaults.numrandomization = 5000;
defaults.minnbchan = 0;
defaults.randomseed = 20260717;
defaults.interpolate = true;

names = fieldnames(defaults);
for i = 1:numel(names)
    if ~isfield(cfg,names{i}) || isempty(cfg.(names{i}))
        cfg.(names{i}) = defaults.(names{i});
    end
end
cfg.metric = char(lower(string(cfg.metric)));
valid_metric = {'signed_delta','magnitude_delta', ...
    'target_signed','target_magnitude'};
assert(ismember(cfg.metric,valid_metric), ...
    'cfg.metric must be one of: %s.',strjoin(valid_metric,', '));
end


function effect = select_effect(target,iti,delta,metric)
switch metric
    case 'signed_delta'
        % [(betaW-betaL)_target - (betaW-betaL)_ITI]
        effect = delta.winner_minus_loser;
    case 'magnitude_delta'
        % Polarity-invariant WTA strength, target minus ITI.
        effect = (abs(target.beta_winner)-abs(target.beta_loser)) - ...
            (abs(iti.beta_winner)-abs(iti.beta_loser));
    case 'target_signed'
        % Target-window signed WTA versus zero.
        effect = target.winner_minus_loser;
    case 'target_magnitude'
        % Target-window polarity-invariant WTA strength versus zero.
        effect = abs(target.beta_winner)-abs(target.beta_loser);
    otherwise
        error('Unknown metric: %s',metric);
end
end


function source = standardize_source(source,grid,name)
assert(isstruct(source) && isfield(source,'pos') && ...
    isfield(source,'inside'),'%s is not a valid source structure.',name);
source = ft_convert_units(source,'mm');
assert(size(source.pos,1)==size(grid.pos,1), ...
    '%s grid count differs from template.',name);
max_error = max(abs(double(source.pos(:))-double(grid.pos(:))));
assert(max_error<1e-5,'%s position mismatch: max error %.6g mm.', ...
    name,max_error);
source.pos = grid.pos;
source.dim = grid.dim;
source.unit = 'mm';
end


function inside = inside_to_logical(value,n_grid)
if islogical(value)
    inside = value(:);
else
    inside = false(n_grid,1);
    inside(double(value(:))) = true;
end
assert(numel(inside)==n_grid,'Inside mask has wrong grid size.');
end


function source = make_ft_source(grid,inside,pow)
source = struct();
source.pos = grid.pos;
source.dim = grid.dim;
source.unit = grid.unit;
if isfield(grid,'coordsys'),source.coordsys = grid.coordsys;end
if isfield(grid,'transform'),source.transform = grid.transform;end
source.inside = logical(inside(:));
source.avg = struct();
source.avg.pow = double(pow(:));
source.dimord = 'pos';
end


function summary = summarize_F_clusters(stat_res,mean_effect,pos,alpha)
summary = table();
if ~isfield(stat_res,'posclusters') || isempty(stat_res.posclusters) || ...
        ~isfield(stat_res,'posclusterslabelmat')
    fprintf('No positive F clusters were formed.\n');
    return;
end

prob = [stat_res.posclusters.prob]';
sig_id = find(prob<=alpha);
if isempty(sig_id)
    fprintf('No cluster survived cluster-level alpha %.4f.\n',alpha);
    return;
end

n = numel(sig_id);
cluster_id = nan(n,1);
cluster_p = nan(n,1);
n_voxel = nan(n,1);
peak_F = nan(n,1);
peak_x = nan(n,1);
peak_y = nan(n,1);
peak_z = nan(n,1);
mean_cluster_effect = nan(n,1);
direction = strings(n,1);

for k = 1:n
    id = sig_id(k);
    idx = find(stat_res.posclusterslabelmat(:)==id);
    [peak_F(k),j] = max(stat_res.stat(idx));
    peak_idx = idx(j);
    cluster_id(k) = id;
    cluster_p(k) = prob(id);
    n_voxel(k) = numel(idx);
    peak_x(k) = pos(peak_idx,1);
    peak_y(k) = pos(peak_idx,2);
    peak_z(k) = pos(peak_idx,3);
    mean_cluster_effect(k) = mean(mean_effect(idx),'omitnan');
    if mean_cluster_effect(k)>0
        direction(k) = "positive";
    elseif mean_cluster_effect(k)<0
        direction(k) = "negative";
    else
        direction(k) = "zero";
    end
end

summary = table(cluster_id,cluster_p,n_voxel,peak_F, ...
    peak_x,peak_y,peak_z,mean_cluster_effect,direction);
summary = sortrows(summary,'cluster_p','ascend');
end


%%
function ExportNIfTI(out_file_base,stat_res)
% =========================================================================
% Export the NIfTI statistic for the first significant cluster after cluster correction
% =========================================================================

% -----------------------------
% basic parameter
% -----------------------------
alpha = 0.05;
cluster_type = 'pos';   % 'pos' or 'neg'

% out_file_base = ['F:\xianliang\exp_data\OT_N2_MEG\MEG_source_408IndependDay\beta_plus\' ...
%     'SigmaV_betacov_F(T-0.50_-0.30)_T(-0.50_-0.30)_FI(0.00_0.20)_TI(0.00_0.20)_u1\all_cluster1'];

ft_path = 'D:\software\matlab_toolbox\fieldtrip-20251218';
template_mri_file = fullfile(ft_path, 'template', 'anatomy', 'single_subj_T1.nii');

% -----------------------------
% 1. Select positive/negative cluster
% -----------------------------
switch lower(cluster_type)
    case 'pos'
        cluster_prob = [stat_res.posclusters.prob];
        cluster_label = stat_res.posclusterslabelmat;

    case 'neg'
        cluster_prob = [stat_res.negclusters.prob];
        cluster_label = stat_res.negclusterslabelmat;

    otherwise
        error('cluster_type 必须是 pos 或 neg。');
end

% -----------------------------
% 2. Find cluster Significant after correction cluster
% -----------------------------
sig_cluster_ids = find(cluster_prob <= alpha);

if isempty(sig_cluster_ids)
    error('没有找到 cluster 校正后 p <= %.3f 的显著 cluster。', alpha);
end

% First significant cluster
target_cluster_id = sig_cluster_ids(1);
target_cluster_p = cluster_prob(target_cluster_id);

fprintf('第一个 cluster 校正显著 cluster: ID = %d, p = %.6f\n', ...
    target_cluster_id, target_cluster_p);

% -----------------------------
% 3. On the original source grid generate cluster mask
% -----------------------------
cluster_mask = cluster_label == target_cluster_id;

if ~any(cluster_mask(:))
    error('labelmat 中没有找到 cluster ID = %d 的点。', target_cluster_id);
end
% ----------------------------
% 4. Retain only this cluster statistical values within
% -----------------------------
source_cluster = [];
% Retain spatial information
source_cluster.pos = stat_res.pos;
source_cluster.dim = stat_res.dim;
source_cluster.inside = stat_res.inside;

if isfield(stat_res, 'unit');source_cluster.unit = stat_res.unit;end
if isfield(stat_res, 'coordsys');source_cluster.coordsys = stat_res.coordsys;end
if isfield(stat_res, 'transform');source_cluster.transform = stat_res.transform;end

% Statistical values
stat_vec = stat_res.stat(:);
mask_vec = cluster_mask(:);

cluster_stat = zeros(size(stat_vec));
cluster_stat(mask_vec) = stat_vec(mask_vec);
cluster_stat(~isfinite(cluster_stat)) = 0;

source_cluster.stat = reshape(cluster_stat, size(stat_res.stat));
source_cluster.statdimord = 'pos';

% -----------------------------
% 5. Interpolate to standard MRI space
% -----------------------------
template_mri = ft_read_mri(template_mri_file);
template_mri = ft_convert_units(template_mri, 'mm');

cfg_interp = [];
cfg_interp.parameter = 'stat';

% cluster label / mask This type of result generally uses nearest.
% Although statistical values are saved here, they are retained only within cluster the region, Therefore nearest better preserves cluster boundary.
cfg_interp.interpmethod = 'nearest';
stat_cluster_interp = ft_sourceinterpolate(cfg_interp, source_cluster, template_mri);
stat_cluster_interp.stat(~isfinite(stat_cluster_interp.stat)) = 0;
% -----------------------------
% 6. Write NIfTI
% -----------------------------
cfg_write = [];
cfg_write.filename = out_file_base;  % No need to add .nii
cfg_write.filetype = 'nifti';
cfg_write.parameter = 'stat';
cfg_write.datatype = 'single';
cfg_write.scaling = 'no';
ft_volumewrite(cfg_write, stat_cluster_interp);
fprintf('✅ 已输出第一个 cluster 校正显著 cluster 的 NIfTI:\n%s.nii\n', out_file_base);
end
