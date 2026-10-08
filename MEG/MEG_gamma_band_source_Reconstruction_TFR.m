%% Grid-wise raw source projection -> source-signal SVD -> ROI TFR
% This script uses the current whole-brain rank-3 DICS definition.
% It does not perform SVD on the ROI filter bank.
%
% Main steps:
%   (1) define a significant functional/OXTR ROI;
%   (2) build one common OT/PL DICS filter at every ROI grid;
%   (3) project short-window raw trials at every grid;
%   (4) estimate one source-space PC1 from the grid-wise source signals;
%   (5) project long-window raw trials and apply the fixed PC1;
%   (6) calculate trial-wise ROI TFR.

clear;

%% Configuration
cfg = struct();
cfg.n_subject = 19;
cfg.filter_time = [-0.5 -0.3];
cfg.tfr_time_high = [-1.2 0.8];
cfg.tfr_time_low  = [-1.8 1.3];
% Final time range retained in the saved TFR structures.
cfg.output_tfr_time = [-1 0.5];
cfg.high_fields_contain_merged_low_high = true;
cfg.tfr_time = [min(cfg.tfr_time_high(1),cfg.tfr_time_low(1)), ...
                max(cfg.tfr_time_high(2),cfg.tfr_time_low(2))];
% Assigned below after the frequency-dependent low-frequency windows are
% defined. The raw trials must cover the full convolution support.
cfg.long_time = [];
cfg.expected_fsample = 250;
cfg.target_foi = [30 100];
cfg.lambda = 0.10;
cfg.gain_pinv_rtol = 1e-10;

% Match the uploaded whole-brain DICS source code.
cfg.apply_day_unit_gain = false;
cfg.grid_ung_for_source_svd = true;
cfg.final_scalar_ung = true;
cfg.tfr_high_freq = 30:2:120;
cfg.tfr_high_win = 0.25*ones(size(cfg.tfr_high_freq));
% The uploaded workflow uses a frequency-dependent window for 2--30 Hz.
cfg.tfr_low_freq = 2:1:30;
cfg.tfr_low_cycles = linspace(3,6,numel(cfg.tfr_low_freq));
cfg.tfr_low_win = cfg.tfr_low_cycles ./ cfg.tfr_low_freq;
cfg.tfr_step = 0.05;
cfg.tfr_trim_percent = 4;
cfg.tfr_edge_padding = 0.5*max([cfg.tfr_high_win cfg.tfr_low_win]) + 0.02;
cfg.long_time = [cfg.tfr_time(1)-cfg.tfr_edge_padding, ...
                 cfg.tfr_time(2)+cfg.tfr_edge_padding];
cfg.tfr_high_smooth = min(10*ones(size(cfg.tfr_high_freq)), ...
    cfg.expected_fsample/2-1-cfg.tfr_high_freq);
assert(all(cfg.tfr_high_smooth>0),'High-frequency TFR smoothing reaches Nyquist.');

cfg.ft_path = 'D:\software\matlab_toolbox\fieldtrip-20251218';
addpath(cfg.ft_path);
ft_defaults;

cfg.template_file = fullfile(cfg.ft_path,'template','sourcemodel', ...
    'standard_sourcemodel3d5mm.mat');
T = load(cfg.template_file,'sourcemodel');
template_grid = ft_convert_units(T.sourcemodel,'mm');

cfg.anatomy_dir = ...
    'F:\xianliang\exp_data\OT_N2_MEG\MRI_Preprocess\head_model';
cfg.data_root = ...
    'F:\xianliang\exp_data\OT_N2_MEG\MEG-Preproocess2\pressPutton';
cfg.lf_root = ...
    'F:\xianliang\exp_data\OT_N2_MEG\MRI_Preprocess\lead_field_pressPutton_3run(250hz)';
cfg.behavior_file = ...
    'E:\xianliang\matlab_m\social_decision_m\data\OT-N2-MEG\MEG_data\pressPutton102(250hz)\data_vdx.mat';
cfg.head_file = ...
    'F:\xianliang\exp_data\OT_N2_MEG\MRI_Preprocess\out_run_4condition\checkHead_idx_2day.mat';

cfg.condition_dirs = {
    fullfile(cfg.data_root,'social_OT_306(250hz)')
    fullfile(cfg.data_root,'social_PL_306(250hz)')
    fullfile(cfg.data_root,'nonsocial_OT_306(250hz)')
    fullfile(cfg.data_root,'nonsocial_PL_306(250hz)')};

cfg.lf_dirs = {
    fullfile(cfg.lf_root,'social-OT(250hz)')
    fullfile(cfg.lf_root,'social-PL(250hz)')
    fullfile(cfg.lf_root,'nonsocial-OT(250hz)')
    fullfile(cfg.lf_root,'nonsocial-PL(250hz)')};

cfg.stat_file = ['F:\xianliang\exp_data\OT_N2_MEG\MEG_source_408IndependDay\' ...
    'unitGain\Source_408DICS_rank3_noDayUnitGain_unitWhiteNoiseGain_30-100Hz_' ...
    '(-0.5_-0.3)\stat_cluster.mat'];
cfg.oxtr_file = ...
    'E:\xianliang\matlab_m\social_decision_m\data\OT-N2-MEG\MEG_data\MEG_mask\article\oxtr_mask\oxtr_50.nii';
cfg.cluster_id = 1;
cfg.oxtr_threshold = 0.5;

cfg.output_dir = fullfile( ...
    ['F:\xianliang\exp_data\OT_N2_MEG\MEG_source_408IndependDay\unitGain\' ...
    'Source_408DICS_rank3_noDayUnitGain_unitWhiteNoiseGain_30-100Hz_(-0.5_-0.3)\source_TFR_single SVD'],'gridwise_raw_sourceSignalSVD');
if ~exist(cfg.output_dir,'dir'), mkdir(cfg.output_dir); end

%% Group metadata and ROI
B = load(cfg.behavior_file,'data_vdx');
data_vdx = B.data_vdx;
H = load(cfg.head_file,'optimal_idx','isolated_bad_run');
mri_dir = dir(fullfile(cfg.anatomy_dir,'*_anatomy.mat'));
assert(numel(mri_dir)>=cfg.n_subject,'Not enough anatomy files.');

S = load(cfg.stat_file,'stat_res');
stat_res = S.stat_res;
oxtr = ft_read_mri(cfg.oxtr_file);
ci = [];
ci.parameter = 'anatomy';
ci.interpmethod = 'nearest';
oxtr_grid = ft_sourceinterpolate(ci,oxtr,template_grid);

FOR GAMMA-BAND SOURCE ∩ OXYTOCIN RECIPTOR EXPRESSION MAP
roi_idx = find(stat_res.posclusterslabelmat(:)==cfg.cluster_id & ...
    oxtr_grid.anatomy(:)>cfg.oxtr_threshold & template_grid.inside(:));
assert(~isempty(roi_idx),'Functional cluster and OXTR mask have no overlap.');
fprintf('ROI grids: %d\n',numel(roi_idx));
% % FOR GAMMA-BAND SOURCE
% roi_idx = find(stat_res.posclusterslabelmat(:)==cfg.cluster_id &  template_grid.inside(:));
% assert(~isempty(roi_idx),'Functional cluster and OXTR mask have no overlap.');
% fprintf('ROI grids: %d\n',numel(roi_idx));


roi_conf = cell(cfg.n_subject,4);
roi_cong = cell(cfg.n_subject,4);
% roi_conf/roi_cong remain the high-frequency trial-level TFR for backward
% compatibility.  The following fields contain both frequency branches.
roi_conf_low = cell(cfg.n_subject,4);
roi_cong_low = cell(cfg.n_subject,4);
high_mean_conf = cell(cfg.n_subject,4);
high_mean_cong = cell(cfg.n_subject,4);
low_mean_conf = cell(cfg.n_subject,4);
low_mean_cong = cell(cfg.n_subject,4);
% All-condition results: concatenate all trials from the four conditions
% first, then apply the same robust trimming and averaging procedure.
high_all_conf = cell(cfg.n_subject,1);
high_all_cong = cell(cfg.n_subject,1);
low_all_conf = cell(cfg.n_subject,1);
low_all_cong = cell(cfg.n_subject,1);
high_all_mean_conf = cell(cfg.n_subject,1);
high_all_mean_cong = cell(cfg.n_subject,1);
low_all_mean_conf = cell(cfg.n_subject,1);
low_all_mean_cong = cell(cfg.n_subject,1);
qc_all = cell(cfg.n_subject,1);

%% Subject loop
for s = 1:cfg.n_subject
    sub_name = regexprep(mri_dir(s).name,'_anatomy.mat','');
    fprintf('\nSubject %d/%d: %s\n',s,cfg.n_subject,sub_name);
    sub_idx = find(strcmp(data_vdx(:,1),sub_name));
    assert(isscalar(sub_idx),'Subject ID is not unique.');

    [D,L] = load_runs_local(sub_name,sub_idx,data_vdx,cfg);
    D_ot = D([1:3 7:9]);
    D_pl = D([4:6 10:12]);
    LOT = L([1:3 7:9]);
    LPL = L([4:6 10:12]);
    run_ot_original = (1:6)';
    run_pl_original = (1:6)';

    bad_ot = finite_integer_vector(H.isolated_bad_run.ot(s));
    bad_pl = finite_integer_vector(H.isolated_bad_run.pl(s));
    D_ot(bad_ot) = [];
    D_pl(bad_pl) = [];
    LOT(bad_ot) = [];
    LPL(bad_pl) = [];

    run_ot_original(bad_ot) = [];
    run_pl_original(bad_pl) = [];
    keep_ot = ~cellfun(@isempty,D_ot);
    keep_pl = ~cellfun(@isempty,D_pl);
    D_ot = D_ot(keep_ot);
    D_pl = D_pl(keep_pl);
    LOT = LOT(keep_ot);
    LPL = LPL(keep_pl);
    run_ot_original = run_ot_original(keep_ot);
    run_pl_original = run_pl_original(keep_pl);

    OT = append_runs_local(D_ot);
    PL = append_runs_local(D_pl);

    iot = find(run_ot_original==double(H.optimal_idx.ot(s)),1);
    ipl = find(run_pl_original==double(H.optimal_idx.pl(s)),1);
    assert(~isempty(iot)&&~isempty(ipl), ...
        'Optimal leadfield points to a removed or empty run.');
    lf_ot = ft_convert_units(LOT{iot},'mm');
    lf_pl = ft_convert_units(LPL{ipl},'mm');

    cells_long = split_cells_local(OT,PL);
    cells_short = cell(4,2);
    for c = 1:4
        for k = 1:2
            cells_short{c,k} = ft_redefinetrial( ...
                struct('toilim',cfg.filter_time),cells_long{c,k});
        end

    end

    [Cot_conf,Fot_conf] = pooled_csd_local( ...
        cells_short([1 3],:),true,cfg);
    [Cot_cong,~] = pooled_csd_local( ...
        cells_short([1 3],:),false,cfg);
    [Cpl_conf,Fpl_conf] = pooled_csd_local( ...
        cells_short([2 4],:),true,cfg);
    [Cpl_cong,~] = pooled_csd_local( ...
        cells_short([2 4],:),false,cfg);

    % The filter covariance pools both task classes and both contexts within
    % each acquisition day.  It is therefore common to all four conditions.
    Cjoint = blkdiag((Cot_conf+Cot_cong)/2, ...
                     (Cpl_conf+Cpl_cong)/2);
    [Wot,Wpl,qot,qpl,Qf] = build_grid_filter_local( ...
        Cjoint,lf_ot,lf_pl,OT.label,PL.label,roi_idx,cfg);

    [u,Qsvd] = train_source_svd_local(cells_short,Wot,Wpl,qot,qpl,cfg);

    for c = 1:4
        if ismember(c,[1 3]), W=Wot; q=qot;
        else, W=Wpl; q=qpl; end
        roi_conf{s,c} = reconstruct_tfr_local( ...
            cells_long{c,1},W,q,u,cfg,'conflict','high');
        roi_cong{s,c} = reconstruct_tfr_local( ...
            cells_long{c,2},W,q,u,cfg,'congruent','high');
        roi_conf_low{s,c} = reconstruct_tfr_local( ...
            cells_long{c,1},W,q,u,cfg,'conflict','low');
        roi_cong_low{s,c} = reconstruct_tfr_local( ...
            cells_long{c,2},W,q,u,cfg,'congruent','low');

        high_mean_conf{s,c} = trim_mean_tfr_local(roi_conf{s,c},cfg.tfr_trim_percent);
        high_mean_cong{s,c} = trim_mean_tfr_local(roi_cong{s,c},cfg.tfr_trim_percent);
        low_mean_conf{s,c} = trim_mean_tfr_local(roi_conf_low{s,c},cfg.tfr_trim_percent);
        low_mean_cong{s,c} = trim_mean_tfr_local(roi_cong_low{s,c},cfg.tfr_trim_percent);
    end

    % Combine all conflict and all congruent trials across the four
    % conditions before trimming. This is a trial-weighted grand condition
    % result, rather than an unweighted mean of four condition means.
    high_all_conf{s} = concat_tfr_trials_local(roi_conf(s,:));
    high_all_cong{s} = concat_tfr_trials_local(roi_cong(s,:));
    low_all_conf{s} = concat_tfr_trials_local(roi_conf_low(s,:));
    low_all_cong{s} = concat_tfr_trials_local(roi_cong_low(s,:));

    high_all_mean_conf{s} = trim_mean_tfr_local( ...
        high_all_conf{s},cfg.tfr_trim_percent);
    high_all_mean_cong{s} = trim_mean_tfr_local( ...
        high_all_cong{s},cfg.tfr_trim_percent);
    low_all_mean_conf{s} = trim_mean_tfr_local( ...
        low_all_conf{s},cfg.tfr_trim_percent);
    low_all_mean_cong{s} = trim_mean_tfr_local( ...
        low_all_cong{s},cfg.tfr_trim_percent);

    Qsvd.subject = sub_name;
    Qsvd.filter = Qf;
    Qsvd.filter_frequency = Fot_conf;
    Qsvd.filter_frequency_pl = Fpl_conf;
    % Qf.grid_idx contains only grids for which both OT and PL filters were
    % successfully constructed and passed the gain checks.
    Qsvd.grid_idx = Qf.grid_idx;
    Qsvd.apply_day_unit_gain = cfg.apply_day_unit_gain;
    Qsvd.grid_ung_for_source_svd = cfg.grid_ung_for_source_svd;
    Qsvd.final_scalar_ung = cfg.final_scalar_ung;
    qc_all{s} = Qsvd;
    fprintf('PC1 explained variance: %.2f%%\n',100*Qsvd.explained);
end

cfg_saved = cfg;
% Condition aliases follow the uploaded naming convention:
% 1=OT-social, 2=PL-social, 3=OT-nonsocial, 4=PL-nonsocial.
% These aliases intentionally point to the trimmed condition means, matching
% the uploaded script's high_* / low_* save convention.  The unaveraged
% trial-level data remain available in roi_conf/roi_cong and their low fields.
high_os_tfr_conf = high_mean_conf(:,1); high_os_tfr_cong = high_mean_cong(:,1);
high_ps_tfr_conf = high_mean_conf(:,2); high_ps_tfr_cong = high_mean_cong(:,2);
high_on_tfr_conf = high_mean_conf(:,3); high_on_tfr_cong = high_mean_cong(:,3);
high_pn_tfr_conf = high_mean_conf(:,4); high_pn_tfr_cong = high_mean_cong(:,4);
low_os_tfr_conf = low_mean_conf(:,1); low_os_tfr_cong = low_mean_cong(:,1);
low_ps_tfr_conf = low_mean_conf(:,2); low_ps_tfr_cong = low_mean_cong(:,2);
low_on_tfr_conf = low_mean_conf(:,3); low_on_tfr_cong = low_mean_cong(:,3);
low_pn_tfr_conf = low_mean_conf(:,4); low_pn_tfr_cong = low_mean_cong(:,4);
high_all_tfr_conf = high_all_mean_conf;
high_all_tfr_cong = high_all_mean_cong;
low_all_tfr_conf = low_all_mean_conf;
low_all_tfr_cong = low_all_mean_cong;

% Merge only the trial-averaged TFRs. Preserve high_mean_* and low_mean_*
% as the original band-specific results; high_* aliases become broadband.
high_alias_conf={high_os_tfr_conf,high_ps_tfr_conf,high_on_tfr_conf,high_pn_tfr_conf};
high_alias_cong={high_os_tfr_cong,high_ps_tfr_cong,high_on_tfr_cong,high_pn_tfr_cong};
low_alias_conf={low_os_tfr_conf,low_ps_tfr_conf,low_on_tfr_conf,low_pn_tfr_conf};
low_alias_cong={low_os_tfr_cong,low_ps_tfr_cong,low_on_tfr_cong,low_pn_tfr_cong};
for c=1:4
    for j=1:numel(high_alias_conf{c})
        high_alias_conf{c}{j}=merge_low_high_tfr_crop_local( ...
            low_alias_conf{c}{j},high_alias_conf{c}{j},cfg.output_tfr_time);
        high_alias_cong{c}{j}=merge_low_high_tfr_crop_local( ...
            low_alias_cong{c}{j},high_alias_cong{c}{j},cfg.output_tfr_time);
    end
end
os_tfr_conf=high_alias_conf{1}; os_tfr_cong=high_alias_cong{1};
ps_tfr_conf=high_alias_conf{2}; ps_tfr_cong=high_alias_cong{2};
on_tfr_conf=high_alias_conf{3}; on_tfr_cong=high_alias_cong{3};
pn_tfr_conf=high_alias_conf{4}; pn_tfr_cong=high_alias_cong{4};
all_tfr_conf = cell(cfg.n_subject,1);
all_tfr_cong = cell(cfg.n_subject,1);

for j = 1:cfg.n_subject
    all_tfr_conf{j} = merge_low_high_tfr_crop_local( ...
        low_all_tfr_conf{j},high_all_tfr_conf{j},cfg.output_tfr_time);

    all_tfr_cong{j} = merge_low_high_tfr_crop_local( ...
        low_all_tfr_cong{j},high_all_tfr_cong{j},cfg.output_tfr_time);
end

subject_ids = cellfun(@(x)x.subject,qc_all,'UniformOutput',false);
save(fullfile(cfg.output_dir,'gridwise_raw_sourceSignalSVD_TFR.mat'), ...
    'os_tfr_conf','os_tfr_cong', ...
    'ps_tfr_conf','ps_tfr_cong', ...
    'on_tfr_conf','on_tfr_cong', ...
    'pn_tfr_conf','pn_tfr_cong', ...
    'all_tfr_conf','all_tfr_cong', ...
    'subject_ids','cfg_saved','roi_idx','-v7.3');

%% ----------------------------- Functions -------------------------------
function [D,L] = load_runs_local(name,si,b,cfg)
D=cell(1,12); L=cell(1,12); n=1;
for c=1:4
    bh=b{si,c+1};
    M=load_one_main_local(dir(fullfile(cfg.condition_dirs{c},[name '*'])));
    assert(numel(M.trial)==numel(bh.is_valid),'MEG/behavior trial mismatch.');
    LF=load(fullfile(dir(fullfile(cfg.lf_dirs{c},[name '*'])).folder, ...
        dir(fullfile(cfg.lf_dirs{c},[name '*'])).name));
    for r=1:3
        ii=find(double(bh.is_valid(:))==1 & double(bh.run_idx(:))==r);
        if isempty(ii),n=n+1;continue;end
        X=ft_selectdata(struct('trials',ii),M);
        X=ft_redefinetrial(struct('toilim',cfg.long_time),X);
        X.trialinfo=[double(bh.is_conflict(ii)),repmat(c,numel(ii),1)];
        D{n}=X; L{n}=LF.leadfield_individual{r}; n=n+1;
    end
end
end

function M=load_one_main_local(d)
assert(isscalar(d),'Expected one MEG file.');
S=load(fullfile(d.folder,d.name)); f=fieldnames(S);
assert(isscalar(f),'MEG MAT must contain one main variable.');
M=S.(f{1});
end

function out=append_runs_local(x)
x=x(~cellfun(@isempty,x)); assert(~isempty(x),'No valid runs.');
% Raw files contain 306 MEG channels plus three auxiliary channels.  The
% leadfields and DICS filters in this pipeline use the 204 planar/gradiometer
% channels only, so remove non-MEG and magnetometer channels before comparing
% labels or appending runs.
for k=1:numel(x)
    cfg_grad=[];
    cfg_grad.channel='MEGGRAD';
    x{k}=ft_selectdata(cfg_grad,x{k});
end

ref=cellstr(x{1}.label(:));
assert(numel(unique(ref))==numel(ref), ...
    'Reference run contains duplicate channel labels.');
for k=1:numel(x)
    lab=cellstr(x{k}.label(:));
    assert(numel(unique(lab))==numel(lab), ...
        'Run %d contains duplicate channel labels.',k);
    [tf,ord]=ismember(ref,lab);
    assert(all(tf) && numel(lab)==numel(ref), ...
        'Run %d has missing or extra channels relative to the reference run.',k);
    if ~isequal(ord(:),(1:numel(ref))')
        % Reorder the numeric trial matrices explicitly.  This avoids a
        % FieldTrip-version-dependent behavior in ft_selectdata, which may
        % select the requested channels but retain their original order.
        for it=1:numel(x{k}.trial)
            X=x{k}.trial{it};
            assert(size(X,1)==numel(lab), ...
                'Run %d trial %d has an invalid channel dimension.',k,it);
            x{k}.trial{it}=X(ord,:);
        end
        x{k}.label=x{k}.label(ord);

        % Keep sensor metadata in the same order when it is present.
        if isfield(x{k},'grad') && isfield(x{k}.grad,'label')
            glab=cellstr(x{k}.grad.label(:));
            [gtf,gord]=ismember(ref,glab);
            if all(gtf) && numel(glab)==numel(ref)
                x{k}.grad.label=x{k}.grad.label(gord);
                x{k}.grad=reorder_grad_fields_local(x{k}.grad,gord);
            end
        end

        % Header channel metadata, when present, follows the data-channel
        % order. Trial/sample metadata are deliberately left unchanged.
        if isfield(x{k},'hdr') && isstruct(x{k}.hdr) && ...
                isfield(x{k}.hdr,'label')
            hlabel=cellstr(x{k}.hdr.label(:));
            [htf,hord]=ismember(ref,hlabel);
            if all(htf) && numel(hlabel)==numel(ref)
                x{k}.hdr.label=x{k}.hdr.label(hord);
                for hf={'chantype','chanunit'}
                    fn=hf{1};
                    if isfield(x{k}.hdr,fn) && ...
                            numel(x{k}.hdr.(fn))==numel(hord)
                        x{k}.hdr.(fn)=x{k}.hdr.(fn)(hord);
                    end
                end
            end
        end
    end
    assert(isequal(cellstr(x{k}.label(:)),ref), ...
        'Explicit channel reordering failed for run %d.',k);
end
out=ft_appenddata(struct('keepsampleinfo','no'),x{:});
end

function grad=reorder_grad_fields_local(grad,ord)
% Reorder row-wise sensor metadata without changing coil geometry itself.
% chanpos/chanori/chantype/chanunit are channel-level fields.
% coilpos and coilori are coil-level fields and must not be reordered by
% channel index merely because their row count happens to match.
row_fields={'chanpos','chanori','chantype','chanunit'};
for i=1:numel(row_fields)
    f=row_fields{i};
    if isfield(grad,f)
        v=grad.(f);
        if size(v,1)==numel(ord)
            grad.(f)=v(ord,:);
        elseif isvector(v) && numel(v)==numel(ord)
            grad.(f)=v(ord);
        end
    end
end
if isfield(grad,'tra') && size(grad.tra,1)==numel(ord)
    grad.tra=grad.tra(ord,:);
end
end

function cells=split_cells_local(OT,PL)
cells=cell(4,2);
for c=1:4
    if ismember(c,[1 3]),D=OT;else,D=PL;end
    % trialinfo(:,1): 1=conflict, 0=congruent
    % trialinfo(:,2): condition code, 1=OT-social, 2=PL-social,
    %                 3=OT-nonsocial, 4=PL-nonsocial
    idx_condition = D.trialinfo(:,2)==c;

    idx_conflict = find(idx_condition & D.trialinfo(:,1)==1);
    idx_congruent = find(idx_condition & D.trialinfo(:,1)==0);

    assert(~isempty(idx_conflict), ...
        'Condition %d has no valid conflict trials.',c);
    assert(~isempty(idx_congruent), ...
        'Condition %d has no valid congruent trials.',c);

    cells{c,1}=ft_selectdata( ...
        struct('trials',idx_conflict),D);
    cells{c,2}=ft_selectdata( ...
        struct('trials',idx_congruent),D);

    assert(all(cells{c,1}.trialinfo(:,2)==c) && ...
           all(cells{c,2}.trialinfo(:,2)==c), ...
        'Condition split failed for condition %d.',c);
end
end

function [C,F]=pooled_csd_local(cells,is_conflict,cfg)
% cells are [social; nonsocial] for one day. Each row contains
% {conflict, congruent}.  For the selected trial type, concatenate all
% social and nonsocial trials first, then estimate one pooled CSD.  This
% matches the whole-brain DICS implementation with pooled-trial weighting.
if is_conflict
    col=1;
else
    col=2;
end
D_social = cells{1,col};
D_nonsocial = cells{2,col};
assert(~isempty(D_social.trial) && ~isempty(D_nonsocial.trial), ...
    'Both social and nonsocial trial sets are required for pooled CSD.');

% append_runs_local also verifies/reorders the 204 MEGGRAD labels.
D_pool = append_runs_local({D_social,D_nonsocial});

fc=[]; fc.method='mtmfft'; fc.output='fourier'; fc.taper='dpss';
fc.channel='MEGGRAD'; fc.foi=mean(cfg.target_foi);
fc.tapsmofrq=diff(cfg.target_foi)/2; fc.keeptrials='yes';
F=ft_freqanalysis(fc,D_pool);
C=csd_local(F);
end

function C=csd_local(F)
X=F.fourierspctrm(:,:,1);
C=(X'*X)/size(X,1);
C=(C+C')/2;
end

function [Wot,Wpl,qot,qpl,Q]=build_grid_filter_local( ...
    C,L1,L2,lab1,lab2,roi,cfg)
[tf1,i1]=ismember(lab1,L1.label(:));
[tf2,i2]=ismember(lab2,L2.label(:));
assert(all(tf1)&&all(tf2),'CSD/leadfield labels do not match.');
scale=real(trace(C))/size(C,1);
invC=pinv((C+C')/2+cfg.lambda*scale*eye(size(C)));
ng=numel(roi); Wot=cell(ng,1);Wpl=cell(ng,1);
qot=nan(ng,1);qpl=nan(ng,1); Q=struct(); Q.n_valid=0;
Q.gain=cell(ng,2);
for k=1:ng
    g=roi(k); A=double(L1.leadfield{g}(i1,:));
    B=double(L2.leadfield{g}(i2,:));
    if isempty(A)||isempty(B),continue;end
    Lj=[A;B]; D=Lj'*invC*Lj;
    W=pinv(D)*(Lj'*invC);
    W1=W(:,1:numel(lab1)); W2=W(:,numel(lab1)+1:end);
    [W1,Q1]=gain_local(W1,A,cfg.apply_day_unit_gain,cfg.gain_pinv_rtol);
    [W2,Q2]=gain_local(W2,B,cfg.apply_day_unit_gain,cfg.gain_pinv_rtol);
    if any(~isfinite([W1(:);W2(:)])),continue;end
    Wot{k}=W1;Wpl{k}=W2;qot(k)=real(trace(W1*W1'));
    qpl(k)=real(trace(W2*W2'));Q.gain{k,1}=Q1;Q.gain{k,2}=Q2;
    Q.n_valid=Q.n_valid+1;
end
ok=~cellfun(@isempty,Wot)&~cellfun(@isempty,Wpl)&qot>0&qpl>0;
assert(any(ok),'No valid grid filter.');
Wot=Wot(ok);Wpl=Wpl(ok);qot=qot(ok);qpl=qpl(ok);
Q.grid_idx=roi(ok);Q.n_valid=sum(ok);Q.rank=3;
end

function [u,Q]=train_source_svd_local(cells,Wot,Wpl,qot,qpl,cfg)
ng=numel(Wot); r=size(Wot{1},1); nf=ng*r;
R=zeros(nf,nf); trc=zeros(4,2);
for c=1:4
    if ismember(c,[1 3]),W=Wot;q=qot;else,W=Wpl;q=qpl;end
    for k=1:2
        D=cells{c,k}; RR=zeros(nf,nf); nobs=0;
        for t=1:numel(D.trial)
            X=double(D.trial{t}); if any(~isfinite(X(:))),continue;end
            Z=grid_project_local(X,W,q,cfg.grid_ung_for_source_svd);
            RR=RR+Z*Z'; nobs=nobs+size(Z,2);
        end
        assert(nobs>0,'Empty source-SVD cell.');
        RR=RR/nobs; trc(c,k)=real(trace(RR)); R=R+RR/8;
    end
end
R=(R+R')/2; [U,d]=eig(R,'vector'); [d,ord]=sort(real(d),'descend');
assert(d(1)>0,'Invalid source second moment.');
u=real(U(:,ord(1))); [~,a]=max(abs(u)); if u(a)<0,u=-u;end
Q.explained=d(1)/sum(max(d,0));Q.eigenvalues=d;
Q.n_grid=ng;Q.n_orientation=r;Q.cell_trace=trc;
Q.method='gridwise raw source projection followed by source-signal SVD';
end

function Z=grid_project_local(X,W,q,useung)
ng=numel(W); r=size(W{1},1); Z=zeros(ng*r,size(X,2));
for k=1:ng
    w=W{k}; if useung,w=w/sqrt(q(k));end
    Z((k-1)*r+(1:r),:)=w*X;
end
end

function T=reconstruct_tfr_local(D,W,q,u,cfg,tag,band)
V=D;V.label={'ROI_sourceSignalSVD'};V.trial=cell(numel(D.trial),1);
for t=1:numel(D.trial)
    Z=grid_project_local(double(D.trial{t}),W,q,cfg.grid_ung_for_source_svd);
    y=real(u'*Z);
    V.trial{t}=y;
end
if cfg.final_scalar_ung
    w_eff = effective_filter_local(W,q,u,cfg.grid_ung_for_source_svd);
    q_eff = real(w_eff*w_eff');
    assert(isfinite(q_eff)&&q_eff>0,'Invalid final scalar white-noise gain.');
    for t=1:numel(V.trial),V.trial{t}=V.trial{t}/sqrt(q_eff);end
end
fc=[];fc.method='mtmconvol';fc.output='pow';fc.channel='all';
switch lower(band)
    case 'high'
        fc.taper='dpss';
        fc.foi=cfg.tfr_high_freq;
        fc.t_ftimwin=cfg.tfr_high_win;
        fc.tapsmofrq=cfg.tfr_high_smooth;
        fc.toi=cfg.tfr_time_high(1):cfg.tfr_step:cfg.tfr_time_high(2);
    case 'low'
        fc.taper='hanning';
        fc.foi=cfg.tfr_low_freq;
        fc.t_ftimwin=cfg.tfr_low_win;
        fc.toi=cfg.tfr_time_low(1):cfg.tfr_step:cfg.tfr_time_low(2);
    otherwise
        error('Unknown TFR band: %s',band);
end
fc.keeptrials='yes';fc.pad='nextpow2';
T=ft_freqanalysis(fc,V);
% FieldTrip cfg.previous can retain a large processing history. Analysis
% settings are saved once in cfg_saved, rather than in every trial TFR.
if isfield(T,'cfg'), T=rmfield(T,'cfg'); end
if ~isfield(T,'trialinfo'), T.trialinfo=D.trialinfo; end
T=crop_tfr_time_local(T,cfg.output_tfr_time);
T.source_svd_method=tag;
T.tfr_band=band;
end

function T=crop_tfr_time_local(T,timewin)
% Crop the already computed TFR without recomputing the convolution.
assert(numel(timewin)==2 && timewin(2)>timewin(1), ...
    'Invalid output TFR time window.');
idx=T.time>=timewin(1)-1e-10 & T.time<=timewin(2)+1e-10;
assert(any(idx),'Output TFR time window has no sampled time points.');
T.time=T.time(idx);
if isfield(T,'powspctrm')
    if ndims(T.powspctrm)==4
        T.powspctrm=T.powspctrm(:,:,:,idx);
    else
        T.powspctrm=T.powspctrm(:,:,idx);
    end
end
end

function Tm=trim_mean_tfr_local(T,pct)
% Robust point-wise MAD/RMS trial rejection, followed by a trial mean.
assert(isfield(T,'powspctrm'),'TFR has no powspctrm.');
assert(strcmp(T.dimord,'rpt_chan_freq_time'), ...
    'Expected trial-level TFR with dimord rpt_chan_freq_time.');
P=double(T.powspctrm);
P=reshape(P,size(P,1),numel(T.label),numel(T.freq),numel(T.time));
n=size(P,1);
[Ptrim,removed_idx]=trim_by_percent_strict(P,1,pct);
kept_idx=setdiff((1:n)',removed_idx(:));
Tm=T;
if isfield(Tm,'cfg'), Tm=rmfield(Tm,'cfg'); end
Tm.dimord='chan_freq_time';
Pm=mean(Ptrim,1,'omitnan');
Tm.powspctrm=reshape(Pm,[numel(T.label) numel(T.freq) numel(T.time)]);
Tm.trim_percent=pct;
Tm.out_idx=removed_idx;
Tm.removed_trial_idx=removed_idx;
Tm.kept_trial_idx=kept_idx;
Tm.n_trials_before=n;
Tm.n_trials_used=numel(kept_idx);
Tm.actual_trim_percent=100*numel(removed_idx)/n;
Tm.trim_method='pointwise median/MAD robust Z-score aggregated by RMS';
trial_fields={'trialinfo','cumtapcnt','cumsumcnt'};
for j=1:numel(trial_fields)
    if isfield(Tm,trial_fields{j}), Tm=rmfield(Tm,trial_fields{j}); end
end
end

function [trimmed_data,removed_idx]=trim_by_percent_strict(data,dim,percent)
% Same scoring/count rule as the uploaded function for percent>0 and n>1.
% Kept local so the reconstruction remains a single .m file.
validateattributes(percent,{'numeric'},{'scalar','real','finite','>=',0,'<=',100});
validateattributes(dim,{'numeric'},{'scalar','integer','>=',1,'<=',ndims(data)});
sz=size(data); n=sz(dim);
assert(n>=1,'Cannot trim an empty trial dimension.');
if percent==0 || n==1
    trimmed_data=data; removed_idx=[]; return;
end
n_remove=max(1,min(round(n*percent/100),n-1));
perm_order=1:ndims(data);
perm_order(perm_order==dim)=[];
perm_order=[dim perm_order];
data_perm=permute(data,perm_order);
matrix_data=reshape(data_perm,n,[]);
assert(~any(isinf(matrix_data(:))),'TFR contains Inf; inspect the input data.');
voxel_med=median(matrix_data,1,'omitnan');
voxel_mad=median(abs(matrix_data-voxel_med),1,'omitnan');
zero_mad_mask=voxel_mad==0;
if any(zero_mad_mask)
    voxel_std=std(matrix_data(:,zero_mad_mask),1,'omitnan');
    voxel_mad(zero_mad_mask)=voxel_std+eps;
end
Z_matrix=abs(matrix_data-voxel_med)./(1.4826*voxel_mad);
slice_scores=sqrt(mean(Z_matrix.^2,2,'omitnan'));
% Explicitly rank completely missing trials as most abnormal.
slice_scores(~isfinite(slice_scores))=Inf;
[~,sorted_idx]=sort(slice_scores,'descend');
removed_idx=sort(sorted_idx(1:n_remove)');
kept_idx=setdiff(1:n,removed_idx);
new_sz=sz; new_sz(dim)=numel(kept_idx);
trimmed_perm=reshape(matrix_data(kept_idx,:),new_sz(perm_order));
[~,inv_perm]=sort(perm_order);
trimmed_data=permute(trimmed_perm,inv_perm);
fprintf('Robust trial trimming: removed %d/%d trials (%.2f%%).\n', ...
    n_remove,n,100*n_remove/n);
end

function T=concat_tfr_trials_local(Tcell)
% Concatenate keeptrials='yes' TFR structures along the trial dimension.
% The input is one subject's four condition-specific TFRs.
Tcell=Tcell(~cellfun(@isempty,Tcell));
assert(~isempty(Tcell),'No condition-level TFR is available.');
T=Tcell{1};
if isfield(T,'cfg'), T=rmfield(T,'cfg'); end
assert(strcmp(T.dimord,'rpt_chan_freq_time'), ...
    'Expected rpt_chan_freq_time TFR for trial concatenation.');
P=T.powspctrm;
for k=2:numel(Tcell)
    Tk=Tcell{k};
    assert(strcmp(Tk.dimord,'rpt_chan_freq_time'), ...
        'All TFR inputs must have rpt_chan_freq_time dimord.');
    assert(isequal(T.label,Tk.label) && isequal(T.freq,Tk.freq) && ...
        isequal(T.time,Tk.time), ...
        'TFR grids differ across conditions.');
    P=cat(1,P,Tk.powspctrm);
end
T.powspctrm=P;
T.dimord='rpt_chan_freq_time';
T.trialinfo=[];
T.n_trials_total=size(P,1);
end

function [Wc,Q]=gain_local(W,L,apply,rtol)
G=W*L;[U,S,V]=svd(G,'econ');s=real(diag(S));
keep=s>rtol*max(s);
if isempty(s)||~any(keep),Wc=nan(size(W));Q=[];return;end
P=V(:,keep)*diag(1./s(keep))*U(:,keep)';
if apply,Wc=P*W;else,Wc=W;end
Q.applied=apply;Q.supported_rank=sum(keep);
Q.error_before=norm(G-P*G,'fro')/max(norm(P*G,'fro'),realmin);
Q.error_after=norm(Wc*L-P*G,'fro')/max(norm(P*G,'fro'),realmin);
end

function v=finite_integer_vector(x)
if iscell(x), x=cell2mat(x(:)); end
v=double(x(:));v=v(isfinite(v)&v==round(v)&v>=1);
end

function j=adjusted_index(x,bad,n)
j=double(x);j=j-sum(bad<j);assert(j>=1&&j<=n,'Bad adjusted leadfield index.');
end

function tfr_merge = merge_low_high_tfr_crop_local(tfr_low, tfr_high, timewin)
% =========================================================================
% Combine low and high frequencies TFR, and crop to the specified time window
%
% New rule: 
%   freq < 31 Hz  -> Use low-frequency TFR
%   freq > 31 Hz  -> Use high-frequency TFR
%   freq == 31 Hz -> Discard
%
% Input requirements: 
%   tfr_low.powspctrm  : chan × freq × time
%   tfr_high.powspctrm : chan × freq × time
%
% Output: 
%   freq Sorted in ascending order
%   powspctrm frequency dimension and freq order must match exactly
% =========================================================================

assert(isfield(tfr_low, 'powspctrm') && isfield(tfr_high, 'powspctrm'), ...
    '输入必须包含 powspctrm。');

assert(isfield(tfr_low, 'freq') && isfield(tfr_high, 'freq'), ...
    '输入必须包含 freq。');

assert(isfield(tfr_low, 'time') && isfield(tfr_high, 'time'), ...
    '输入必须包含 time。');

assert(numel(timewin) == 2 && timewin(2) > timewin(1), ...
    'timewin 必须是 [start end]，例如 [-1.2 0.8]。');

if isfield(tfr_low, 'dimord')
    assert(strcmp(tfr_low.dimord, 'chan_freq_time'), ...
        '当前函数只支持 tfr_low.dimord = chan_freq_time。');
end

if isfield(tfr_high, 'dimord')
    assert(strcmp(tfr_high.dimord, 'chan_freq_time'), ...
        '当前函数只支持 tfr_high.dimord = chan_freq_time。');
end

assert(isequal(tfr_low.label(:),tfr_high.label(:)), ...
    'Low/high TFR channel labels or order differ.');

tol = 1e-6;
cut_freq = 31;

low_freq_raw  = tfr_low.freq(:)';
high_freq_raw = tfr_high.freq(:)';

% ---------------------------------------------------------
% 0. Select frequencies by the rule: 
%    Retain only low-frequency <30, Retain only high-frequency >30
% ---------------------------------------------------------
low_keep_raw  = low_freq_raw  < cut_freq - tol;
high_keep_raw = high_freq_raw > cut_freq + tol;

if ~any(low_keep_raw)
    error('低频 TFR 中没有 freq < %.1f Hz 的频点。当前 low freq 为：%s', ...
        cut_freq, mat2str(low_freq_raw));
end

if ~any(high_keep_raw)
    error('高频 TFR 中没有 freq > %.1f Hz 的频点。当前 high freq 为：%s', ...
        cut_freq, mat2str(high_freq_raw));
end

low_freq_keep_raw  = low_freq_raw(low_keep_raw);
high_freq_keep_raw = high_freq_raw(high_keep_raw);

% ---------------------------------------------------------
% 1. Check for duplicates within retained frequencies
% ---------------------------------------------------------
low_freq_check = sort(low_freq_keep_raw);
high_freq_check = sort(high_freq_keep_raw);

if any(diff(low_freq_check) <= tol)
    error('保留后的 low freq 内部存在重复频点：%s', ...
        mat2str(low_freq_keep_raw));
end

if any(diff(high_freq_check) <= tol)
    error('保留后的 high freq 内部存在重复频点：%s', ...
        mat2str(high_freq_keep_raw));
end

% ---------------------------------------------------------
% 2. Retain low-frequency <30 Hz, and sort ascending; Reorder concurrently power
% ---------------------------------------------------------
[low_freq, low_sort_idx] = sort(low_freq_keep_raw, 'ascend');
low_keep_idx = find(low_keep_raw);
low_freq_idx_sorted = low_keep_idx(low_sort_idx);

low_pow_sorted = double(tfr_low.powspctrm(:, low_freq_idx_sorted, :));

% ---------------------------------------------------------
% 3. Retain high-frequency >30 Hz, and sort ascending; Reorder concurrently power
% ---------------------------------------------------------
[high_freq, high_sort_idx] = sort(high_freq_keep_raw, 'ascend');
high_keep_idx = find(high_keep_raw);
high_freq_idx_sorted = high_keep_idx(high_sort_idx);

high_time_idx = find(tfr_high.time >= timewin(1)-tol & ...
                     tfr_high.time <= timewin(2)+tol);

if isempty(high_time_idx)
    error('高频 TFR 中没有落在目标时间窗 [%.3f %.3f] 内的时间点。', ...
        timewin(1), timewin(2));
end

high_pow_crop = double(tfr_high.powspctrm(:, high_freq_idx_sorted, high_time_idx));

% ---------------------------------------------------------
% 4. Use the high-frequency time axis within timewin the time interval
% ---------------------------------------------------------
merge_time = tfr_high.time(high_time_idx);
merge_time = merge_time(:)';

% ---------------------------------------------------------
% 5. Low-frequency power Interpolate onto the high-frequency time axis
% ---------------------------------------------------------
low_time = tfr_low.time(:)';

if min(merge_time) < min(low_time)-tol || max(merge_time) > max(low_time)+tol
    error(['低频 TFR 时间范围不能覆盖最终合并时间窗。' ...
           '低频范围为 [%.3f %.3f]，目标范围为 [%.3f %.3f]。'], ...
        min(low_time), max(low_time), min(merge_time), max(merge_time));
end

n_chan = size(low_pow_sorted, 1);
n_low_freq = size(low_pow_sorted, 2);
n_low_time = size(low_pow_sorted, 3);

low_pow_2d = reshape(low_pow_sorted, n_chan * n_low_freq, n_low_time)';
low_pow_interp_2d = interp1(low_time, low_pow_2d, merge_time, 'linear');

if any(isnan(low_pow_interp_2d(:)))
    warning('低频插值后出现 NaN，请检查 low_time 是否完整覆盖 merge_time。');
end

low_pow_crop = reshape(low_pow_interp_2d', ...
    n_chan, n_low_freq, numel(merge_time));

% ---------------------------------------------------------
% 6. Concatenate the frequency dimension
% ---------------------------------------------------------
merge_freq_pre = [low_freq, high_freq];
merge_pow_pre = cat(2, low_pow_crop, high_pow_crop);

% ---------------------------------------------------------
% 7. Sort globally again after concatenation, and reorder concurrently power
% ---------------------------------------------------------
[merge_freq, merge_sort_idx] = sort(merge_freq_pre, 'ascend');
merge_pow = merge_pow_pre(:, merge_sort_idx, :);

% ---------------------------------------------------------
% 8. Final checks: Frequency order, power Dimension, 30 Hz Whether it was removed
% ---------------------------------------------------------
if any(diff(merge_freq) <= tol)
    error('合并后频率轴仍存在重复或倒序。当前 merge_freq 为：%s', ...
        mat2str(merge_freq));
end

if any(abs(merge_freq - cut_freq) <= tol)
    error('合并后仍包含 %.1f Hz，请检查频率筛选逻辑。当前 merge_freq 为：%s', ...
        cut_freq, mat2str(merge_freq));
end

if any(merge_freq < cut_freq - tol) && any(merge_freq > cut_freq + tol)
    % Normal: Both frequency bands contain
else
    warning('合并后似乎只包含低频或只包含高频，请检查输入数据。');
end

if size(merge_pow, 2) ~= numel(merge_freq)
    error('powspctrm 的频率维度数量与 freq 数量不一致。');
end

if any(isnan(merge_pow(:)))
    warning('合并后的 powspctrm 中存在 NaN。');
end

% ---------------------------------------------------------
% 9. Organize output structure
% ---------------------------------------------------------
tfr_merge = [];
tfr_merge.label = tfr_low.label;
tfr_merge.dimord = 'chan_freq_time';
tfr_merge.freq = merge_freq;
tfr_merge.time = merge_time;
tfr_merge.powspctrm = merge_pow;

% Auxiliary information
tfr_merge.low_freq = low_freq;
tfr_merge.high_freq = high_freq;
tfr_merge.cut_freq = cut_freq;
tfr_merge.freq_rule = ...
    'freq < 31 Hz uses low-frequency TFR; freq > 31 Hz uses high-frequency TFR; freq == 31 Hz is removed.';

tfr_merge.original_low_freq = low_freq_raw;
tfr_merge.original_high_freq = high_freq_raw;
tfr_merge.original_low_time = tfr_low.time;
tfr_merge.original_high_time = tfr_high.time;
tfr_merge.merge_timewin = timewin;

tfr_merge.merge_note = ...
    'Merged low+high TFR. Low keeps freq < 31 Hz, high keeps freq > 31 Hz. Frequencies are sorted ascending and powspctrm is reordered with the same index.';

if isfield(tfr_low, 'out_idx')
    tfr_merge.low_out_idx = tfr_low.out_idx;
end

if isfield(tfr_high, 'out_idx')
    tfr_merge.high_out_idx = tfr_high.out_idx;
end

fprintf('✅ 合并完成：低频 %.1f–%.1f Hz，高频 %.1f–%.1f Hz；最终 freq %.1f–%.1f Hz，共 %d 个频点；time %.3f–%.3f s，共 %d 个时间点。\n', ...
    min(low_freq), max(low_freq), ...
    min(high_freq), max(high_freq), ...
    min(tfr_merge.freq), max(tfr_merge.freq), numel(tfr_merge.freq), ...
    min(tfr_merge.time), max(tfr_merge.time), numel(tfr_merge.time));

end

function w=effective_filter_local(W,q,u,useung)
ng=numel(W); r=size(W{1},1); nc=size(W{1},2);
Wflat=zeros(ng*r,nc);
for k=1:ng
    wk=W{k};
    if useung, wk=wk/sqrt(q(k)); end
    Wflat((k-1)*r+(1:r),:)=wk;
end
w=u(:)'*Wflat;
end
