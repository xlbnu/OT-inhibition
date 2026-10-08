% iti_delay_num=[0 0.1 0.2];
% for i_delay=1:length(iti_delay_num)
cfg_main = struct();

% Construct beta covariance/filter task time windows
cfg_main.source_cfg.filter_target_window = [-0.5 -0.3];

% Use the constructed filter projection time range
cfg_main.source_cfg.projection_target_window = [-0.5 -0.3];

% filter Used ITI
cfg_main.source_cfg.filter_iti_delay = 0.0;%iti_delay_num(i_delay);
cfg_main.source_cfg.filter_iti_duration = 0.2;
cfg_main.source_cfg.filter_min_iti_duration = 0.19;

% Used for projection ITI
cfg_main.source_cfg.projection_iti_delay = 0.0;%iti_delay_num(i_delay);
cfg_main.source_cfg.projection_iti_duration = 0.2;
cfg_main.source_cfg.projection_min_iti_duration = 0.19;


% Two-dimensional source output: 
%   'magnitude' = sqrt(q1.^2 + q2.^2)
%   'sumsq'     = q1.^2 + q2.^2
%   'both'      = Main/Compatibility fields use magnitude, and additionally save both explicit results
cfg_main.source_cfg.vector_output_mode = 'both';

cfg_main.output_dir = ...
    fullfile('F:\xianliang\exp_data\OT_N2_MEG\MEG_source_408IndependDay\beta_plus',sprintf('SigmaV_betacov_F(T%0.2f_%0.2f)_T(%0.2f_%0.2f)_FI(%0.2f_%0.2f)_TI(%0.2f_%0.2f)_u1',...
    cfg_main.source_cfg.filter_target_window(1),cfg_main.source_cfg.filter_target_window(2),...
    cfg_main.source_cfg.projection_target_window(1),cfg_main.source_cfg.projection_target_window(2),...
    cfg_main.source_cfg.filter_iti_delay,...
    cfg_main.source_cfg.filter_iti_delay+cfg_main.source_cfg.filter_iti_duration,...
    cfg_main.source_cfg.projection_iti_delay,...
    cfg_main.source_cfg.projection_iti_delay+cfg_main.source_cfg.projection_iti_duration));




%%
cfg_main = set_main_defaults(cfg_main);

addpath(cfg_main.ft_path);
ft_defaults;
if ~exist(cfg_main.output_dir,'dir'), mkdir(cfg_main.output_dir); end

%% Behavior, head movement and template grid
behavior_file = load(cfg_main.behavior_file);
required_sets = {'ots_meg','pls_meg','otn_meg','pln_meg'};
for i = 1:numel(required_sets)
    assert(isfield(behavior_file,required_sets{i}), ...
        'Behavior file lacks %s.',required_sets{i});
end
data_bhv = {behavior_file.ots_meg;behavior_file.pls_meg; ...
    behavior_file.otn_meg;behavior_file.pln_meg};
subject_ids = get_and_verify_subject_ids(data_bhv);
n_subject = numel(subject_ids);
if ~isempty(cfg_main.num_subjects)
    assert(cfg_main.num_subjects<=n_subject);
    n_subject = cfg_main.num_subjects;
    subject_ids = subject_ids(1:n_subject);
end

head_mov = load(cfg_main.head_movement_file, ...
    'optimal_idx','isolated_bad_run');
assert(isfield(head_mov,'optimal_idx') && ...
    isfield(head_mov,'isolated_bad_run'), ...
    'Head-movement file lacks optimal_idx/isolated_bad_run.');
head_mov.isolated_bad_run.ot{18}=[3 6];
source_grid = load_source_grid(cfg_main.template_grid_file);
source_grid = ft_convert_units(source_grid,'mm');


%% Subjects
subject_result = struct();
subject_result.subject_ids = subject_ids(:);
subject_result.effect_definition = [ ...
    'SigmaV = chosen + unchosen; DeltaV = chosen - unchosen; ' ...
    'filter and source maps use beta_SigmaV only'];
subject_result.vector_output_mode = ...
    cfg_main.source_cfg.vector_output_mode;
subject_result.quality = cell(n_subject,1);
for i_sub = 1:n_subject
    subject_id = subject_ids{i_sub};
    fprintf('\n============================================================\n');
    fprintf('[%d/%d] : %s\n',i_sub,n_subject,subject_id);

    [ot_data,pl_data,bhv_ot,bhv_pl,ot_lf,pl_lf,quality] = ...
        prepare_subject(i_sub,subject_id,data_bhv,head_mov,cfg_main);

    cfg_beta = struct();

    % Construct a common filter 's target time window
    cfg_beta.filter_target_window = ...
        cfg_main.source_cfg.filter_target_window;

    % Construct a common filter 's ITI
    cfg_beta.filter_iti_delay = ...
        cfg_main.source_cfg.filter_iti_delay;
    cfg_beta.filter_iti_duration = ...
        cfg_main.source_cfg.filter_iti_duration;
    cfg_beta.filter_min_iti_duration = ...
        cfg_main.source_cfg.filter_min_iti_duration;

    % Finally produce target source map time window
    cfg_beta.projection_target_window = ...
        cfg_main.source_cfg.projection_target_window;

    % Finally produce ITI source map time window
    cfg_beta.projection_iti_delay = ...
        cfg_main.source_cfg.projection_iti_delay;
    cfg_beta.projection_iti_duration = ...
        cfg_main.source_cfg.projection_iti_duration;
    cfg_beta.projection_min_iti_duration = ...
        cfg_main.source_cfg.projection_min_iti_duration;

    % beta covariance
    cfg_beta.beta_covariance_mode = 'secondmoment';
    cfg_beta.beta_smooth_samples = 1;

    % rank-2 joint LCMV
    cfg_beta.source_rank = 2;
    cfg_beta.lambda = 0.10;
    cfg_beta.apply_day_gain_calibration = true;
    cfg_beta.min_rcond = 1e-10;
    cfg_beta.vector_output_mode = ...
        cfg_main.source_cfg.vector_output_mode;

    % trial choice
    cfg_beta.condition_include = 1:4;
    cfg_beta.require_conflict = true;
    cfg_beta.min_trial = 20;

    % With the current sensor GLM remain consistent
    cfg_beta.add_change_nuisance = false;
    cfg_beta.add_condition_nuisance = false;
    cfg_beta.add_run_nuisance = false;

    % Whole-brain stage
    cfg_beta.grid_indices = [];
    cfg_beta.keepfilter = false;
    cfg_beta.filter_storage_class = 'single';

    vector_result = run_sigmaV_timebeta_vector_lcmv_core( ...
        ot_data, pl_data, ...
        bhv_ot, bhv_pl, ...
        ot_lf, pl_lf, ...
        source_grid, cfg_beta);

    % Target/noise maps of the SigmaV beta magnitude.
    subject_result.target{i_sub,1} = vector_result.maps.target_unitnoisegain;
    subject_result.iti{i_sub,1} = vector_result.maps.iti_unitnoisegain;
    subject_result.target_minus_iti{i_sub,1} = vector_result.maps.target_minus_iti_unitnoisegain;

    subject_result.target_raw{i_sub,1} = vector_result.maps.target_raw;
    subject_result.iti_raw{i_sub,1} = vector_result.maps.iti_raw;
    subject_result.target_minus_iti_raw{i_sub,1} = vector_result.maps.target_minus_iti_raw;

    % 'both'In this mode, Additionally retain explicitly named magnitude and sumsq results.
    if strcmp(vector_result.cfg.vector_output_mode,'both')
        subject_result.by_measure.magnitude.target_unitnoisegain{i_sub,1} = ...
            vector_result.maps_magnitude.target_unitnoisegain;
        subject_result.by_measure.magnitude.iti_unitnoisegain{i_sub,1} = ...
            vector_result.maps_magnitude.iti_unitnoisegain;
        subject_result.by_measure.magnitude.target_minus_iti_unitnoisegain{i_sub,1} = ...
            vector_result.maps_magnitude.target_minus_iti_unitnoisegain;
        subject_result.by_measure.magnitude.target_raw{i_sub,1} = ...
            vector_result.maps_magnitude.target_raw;
        subject_result.by_measure.magnitude.iti_raw{i_sub,1} = ...
            vector_result.maps_magnitude.iti_raw;
        subject_result.by_measure.magnitude.target_minus_iti_raw{i_sub,1} = ...
            vector_result.maps_magnitude.target_minus_iti_raw;

        subject_result.by_measure.sumsq.target_unitnoisegain{i_sub,1} = ...
            vector_result.maps_sumsq.target_unitnoisegain;
        subject_result.by_measure.sumsq.iti_unitnoisegain{i_sub,1} = ...
            vector_result.maps_sumsq.iti_unitnoisegain;
        subject_result.by_measure.sumsq.target_minus_iti_unitnoisegain{i_sub,1} = ...
            vector_result.maps_sumsq.target_minus_iti_unitnoisegain;
        subject_result.by_measure.sumsq.target_raw{i_sub,1} = ...
            vector_result.maps_sumsq.target_raw;
        subject_result.by_measure.sumsq.iti_raw{i_sub,1} = ...
            vector_result.maps_sumsq.iti_raw;
        subject_result.by_measure.sumsq.target_minus_iti_raw{i_sub,1} = ...
            vector_result.maps_sumsq.target_minus_iti_raw;
    end

    % Separate nGrid x nTime target and ITI SigmaV-beta trajectories.
    % subject_result.dynamic{i_sub,1} = timebeta_result.dynamic;
    subject_result.time{i_sub,1} = vector_result.time;
    subject_result.sensor_beta_time{i_sub,1} = vector_result.sensor_beta_time;
    subject_result.quality{i_sub,1} = struct( ...
        'preparation',quality,'source',vector_result.quality);
    % subject_result.cfg{i_sub,1} = vector_result.cfg;

end

switch cfg_main.source_cfg.vector_output_mode
    case 'magnitude'
        output_measure_tag = 'vectorMagnitude';
    case 'sumsq'
        output_measure_tag = 'vectorSumSq';
    case 'both'
        output_measure_tag = 'vectorMagnitude_and_SumSq';
    otherwise
        error('Unknown vector_output_mode: %s', ...
            cfg_main.source_cfg.vector_output_mode);
end
subject_file = fullfile(cfg_main.output_dir, ...
    sprintf('SigmaV_%s_source.mat',output_measure_tag));

save(subject_file, ...
    'subject_result','cfg_main','-v7.3');

% end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% functions
function [ot_data,pl_data,bhv_ot,bhv_pl,ot_lf,pl_lf,quality] = ...
    prepare_subject(i_sub,subject_id,data_bhv,head_mov,cfg_main)
all_run_data = cell(12,1);
all_run_bhv = cell(12,1);
all_run_meta = cell(12,1);
run_counter = 0;
ot_slots = [1:3 7:9];
pl_slots = [4:6 10:12];
bad_ot = get_subject_index_vector(head_mov.isolated_bad_run.ot,i_sub);
bad_pl = get_subject_index_vector(head_mov.isolated_bad_run.pl,i_sub);
validate_day_run_indices(bad_ot,'OT isolated bad run');
validate_day_run_indices(bad_pl,'PL isolated bad run');
optimal_ot = get_subject_scalar(head_mov.optimal_idx.ot,i_sub);
optimal_pl = get_subject_scalar(head_mov.optimal_idx.pl,i_sub);
validate_day_run_indices(optimal_ot,'OT optimal run');
validate_day_run_indices(optimal_pl,'PL optimal run');
assert(~ismember(optimal_ot,bad_ot), ...
    'OT optimal leadfield run is also marked as an isolated bad run.');
assert(~ismember(optimal_pl,bad_pl), ...
    'PL optimal leadfield run is also marked as an isolated bad run.');
optimal_ot_slot = ot_slots(optimal_ot);
optimal_pl_slot = pl_slots(optimal_pl);
ot_lf = [];
pl_lf = [];

for i_condition = 1:4
    raw_file = find_unique_subject_file( ...
        cfg_main.condition_dirs{i_condition},subject_id);
    data_process = load_raw_with_rejection(raw_file);
    bhv0 = data_bhv{i_condition}(i_sub);
    assert(strcmp(normalize_subject_id(bhv0.subsName),subject_id), ...
        'Behavior subject mismatch in condition %d.',i_condition);
    [bhv0,bad_global,bad_by_run] = remove_rejected_behavior_trials( ...
        bhv0,data_process.rejection_info,cfg_main.trials_per_run);
    n_meg = numel(data_process.trial);
    n_bhv = get_behavior_trial_count(bhv0);
    assert(n_meg==n_bhv, ...
        ['After rejection_info, MEG has %d trials but behavior has %d: ' ...
        '%s condition %d.'],n_meg,n_bhv,subject_id,i_condition);

    bhv0.run_id = make_run_id_after_rejection( ...
        cfg_main.trials_per_run,bad_global);
    bhv0.condition = repmat(i_condition,n_bhv,1);
    bhv0.original_trial_id = setdiff( ...
        (1:3*cfg_main.trials_per_run)',bad_global,'stable');
    bhv0.run_uid = (i_condition-1)*3+bhv0.run_id;

    % rejection_info has already been consumed above. It is not a
    % FieldTrip data dimension, so remove it before FieldTrip functions to
    % prevent repeated "could not determine dimord" warnings.
    data_process = rmfield_if_exists(data_process,{'rejection_info'});

    % cfg_crop = [];
    % cfg_crop.toilim = cfg_main.analysis_window;
    % data_process = ft_redefinetrial(cfg_crop,data_process);

    lf_file = find_unique_subject_file( ...
        cfg_main.leadfield_dirs{i_condition},subject_id);
    condition_slots = (i_condition-1)*3+(1:3);
    needs_lf = ismember(optimal_ot_slot,condition_slots) || ...
        ismember(optimal_pl_slot,condition_slots);
    if needs_lf
        [leadfield_3run,~] = load_three_run_leadfields(lf_file);
    else
        leadfield_3run = [];
    end

    for i_run = 1:3
        run_counter = run_counter+1;
        run_idx = find(bhv0.run_id==i_run & bhv0.conflict==1);
        if run_counter==optimal_ot_slot
            ot_lf = ft_convert_units(leadfield_3run{i_run},'mm');
        elseif run_counter==optimal_pl_slot
            pl_lf = ft_convert_units(leadfield_3run{i_run},'mm');
        end
        all_run_meta{run_counter} = struct( ...
            'condition',i_condition,'run',i_run,'raw_file',raw_file, ...
            'leadfield_file',lf_file,'bad_trials_local',bad_by_run{i_run}, ...
            'n_conflict_trial',numel(run_idx));
        if isempty(run_idx)
            warning('%s condition %d run %d has no retained conflict trial.', ...
                subject_id,i_condition,i_run);
            continue;
        end
        cfg_sel = [];
        cfg_sel.trials = run_idx;
        cfg_sel.channel = 'MEGGRAD';
        run_data = ft_selectdata(cfg_sel,data_process);
        run_data = rmfield_if_exists(run_data, ...
            {'grad','elec','hdr','sampleinfo','rejection_info'});
        run_bhv = behavior_struct_to_table(bhv0,run_idx, ...
            cfg_main.behavior_fields);
        assert(numel(run_data.trial)==height(run_bhv));
        assert(all(run_bhv.conflict==1) && all(run_bhv.run_id==i_run));
        all_run_data{run_counter} = run_data;
        all_run_bhv{run_counter} = run_bhv;
    end
end
assert(run_counter==12);

ot_keep = setdiff(1:6,bad_ot,'stable');
pl_keep = setdiff(1:6,bad_pl,'stable');
[ot_data,bhv_ot,ot_used_slots] = append_day_runs( ...
    all_run_data,all_run_bhv,ot_slots,ot_keep,'OT');
[pl_data,bhv_pl,pl_used_slots] = append_day_runs( ...
    all_run_data,all_run_bhv,pl_slots,pl_keep,'PL');
assert(~isempty(ot_lf) && ~isempty(pl_lf));

quality = struct();
quality.subject_id = subject_id;
quality.optimal_ot_index_within_day = optimal_ot;
quality.optimal_pl_index_within_day = optimal_pl;
quality.optimal_ot_slot = optimal_ot_slot;
quality.optimal_pl_slot = optimal_pl_slot;
quality.optimal_ot_condition = ceil(optimal_ot_slot/3);
quality.optimal_pl_condition = ceil(optimal_pl_slot/3);
quality.optimal_ot_run_within_condition = mod(optimal_ot_slot-1,3)+1;
quality.optimal_pl_run_within_condition = mod(optimal_pl_slot-1,3)+1;
quality.optimal_ot_leadfield_file = ...
    all_run_meta{optimal_ot_slot}.leadfield_file;
quality.optimal_pl_leadfield_file = ...
    all_run_meta{optimal_pl_slot}.leadfield_file;
quality.removed_ot_run_indices = bad_ot;
quality.removed_pl_run_indices = bad_pl;
quality.ot_used_slots = ot_used_slots;
quality.pl_used_slots = pl_used_slots;
quality.run_metadata = all_run_meta;
quality.n_trial_ot = height(bhv_ot);
quality.n_trial_pl = height(bhv_pl);
quality.demean_input_expected = cfg_main.input_is_whole_trial_demeaned;
end


function cfg = set_main_defaults(cfg)
defaults = struct();
defaults.num_subjects = 19;
defaults.ft_path = 'D:\software\matlab_toolbox\fieldtrip-20251218';
defaults.base_data_path = ...
    'F:\xianliang\exp_data\OT_N2_MEG\MEG-Preproocess2\pressPutton';
% defaults.base_data_path = ...
%     'F:\xianliang\exp_data\OT_N2_MEG\MEG-Preproocess2\pressPutton_T(-7_1)s';
defaults.behavior_file = ...
    'E:\xianliang\matlab_m\social_decision_m\data\OT-N2-MEG\behavior\analysis_data\behaviour_meg.mat';
defaults.base_leadfield_path = ...
    'F:\xianliang\exp_data\OT_N2_MEG\MRI_Preprocess\lead_field_pressPutton_3run(250hz)';
defaults.head_movement_file = ...
    'F:\xianliang\exp_data\OT_N2_MEG\MRI_Preprocess\out_run_4condition\checkHead_idx_2day.mat';
defaults.template_grid_file = fullfile(defaults.ft_path,'template', ...
    'sourcemodel','standard_sourcemodel3d5mm.mat');
defaults.output_dir = fullfile( ...
    'F:\xianliang\exp_data\OT_N2_MEG\MEG_source_408IndependDay', ...
    'SigmaV_beta_jointLCMV');
defaults.analysis_window = [-3 3];
defaults.trials_per_run = 80;
defaults.input_is_whole_trial_demeaned = true;
defaults.behavior_fields = {'conflict','ischange','lo_s1','lo_o1', ...
    'condition','run_id','run_uid','oiti','osti2','rtime2', ...
    'original_trial_id','lo_chosen','lo_unchosen','osti'};
names = fieldnames(defaults);
for i = 1:numel(names)
    if ~isfield(cfg,names{i}) || isempty(cfg.(names{i}))
        cfg.(names{i}) = defaults.(names{i});
    end
end

cfg.condition_dirs = {
    fullfile(cfg.base_data_path,'social_OT_306(250hz)')
    fullfile(cfg.base_data_path,'social_PL_306(250hz)')
    fullfile(cfg.base_data_path,'nonsocial_OT_306(250hz)')
    fullfile(cfg.base_data_path,'nonsocial_PL_306(250hz)')};
cfg.leadfield_dirs = {
    fullfile(cfg.base_leadfield_path,'social-OT(250hz)')
    fullfile(cfg.base_leadfield_path,'social-PL(250hz)')
    fullfile(cfg.base_leadfield_path,'nonsocial-OT(250hz)')
    fullfile(cfg.base_leadfield_path,'nonsocial-PL(250hz)')};

if ~isfield(cfg,'source_cfg') || isempty(cfg.source_cfg)
    cfg.source_cfg = struct();
end
source_defaults = struct( ...
    'filter_target_window',[-0.7 -0.3], ...
    'projection_target_window',[-0.5 -0.3], ...
    'filter_iti_duration',0.2, ...
    'filter_iti_delay',0.3, ...
    'filter_min_iti_duration',0.19, ...
    'projection_iti_duration',0.2, ...
    'projection_iti_delay',0.3, ...
    'projection_min_iti_duration',0.19, ...
    'channel_selection','MEGGRAD', ...
    'center_covariance',true, ...
    'lambda',0.10, ...
    'regularization_mode','joint_global', ...
    'apply_day_gain_calibration',true, ...
    'vector_output_mode','magnitude', ...
    'log_transform',false, ...
    'min_trial',20, ...
    'min_rcond',1e-10, ...
    'keepfilter',false, ...
    'filter_storage_class','single');
names = fieldnames(source_defaults);
for i = 1:numel(names)
    if ~isfield(cfg.source_cfg,names{i}) || ...
            isempty(cfg.source_cfg.(names{i}))
        cfg.source_cfg.(names{i}) = source_defaults.(names{i});
    end
end
% Optional independent windows. Old cfg.source_cfg.target_window and
% cfg.source_cfg.iti_* fields remain valid and are used as fallbacks.
if ~isfield(cfg.source_cfg,'filter_target_window') || ...
        isempty(cfg.source_cfg.filter_target_window)
    cfg.source_cfg.filter_target_window = cfg.source_cfg.target_window;
end
if ~isfield(cfg.source_cfg,'projection_target_window') || ...
        isempty(cfg.source_cfg.projection_target_window)
    cfg.source_cfg.projection_target_window = cfg.source_cfg.target_window;
end
if ~isfield(cfg.source_cfg,'filter_iti_duration') || ...
        isempty(cfg.source_cfg.filter_iti_duration)
    cfg.source_cfg.filter_iti_duration = cfg.source_cfg.iti_duration;
end
if ~isfield(cfg.source_cfg,'filter_iti_delay') || ...
        isempty(cfg.source_cfg.filter_iti_delay)
    cfg.source_cfg.filter_iti_delay = cfg.source_cfg.iti_delay;
end
if ~isfield(cfg.source_cfg,'filter_min_iti_duration') || ...
        isempty(cfg.source_cfg.filter_min_iti_duration)
    cfg.source_cfg.filter_min_iti_duration = ...
        cfg.source_cfg.min_iti_duration;
end
if ~isfield(cfg.source_cfg,'projection_iti_duration') || ...
        isempty(cfg.source_cfg.projection_iti_duration)
    cfg.source_cfg.projection_iti_duration = cfg.source_cfg.iti_duration;
end
if ~isfield(cfg.source_cfg,'projection_iti_delay') || ...
        isempty(cfg.source_cfg.projection_iti_delay)
    cfg.source_cfg.projection_iti_delay = cfg.source_cfg.iti_delay;
end
if ~isfield(cfg.source_cfg,'projection_min_iti_duration') || ...
        isempty(cfg.source_cfg.projection_min_iti_duration)
    cfg.source_cfg.projection_min_iti_duration = ...
        cfg.source_cfg.min_iti_duration;
end

end


function ids = get_and_verify_subject_ids(data_bhv)
n = numel(data_bhv{1});
ids = cell(n,1);
for i_sub = 1:n
    ids{i_sub} = normalize_subject_id(data_bhv{1}(i_sub).subsName);
    for i_cond = 2:4
        assert(strcmp(ids{i_sub},normalize_subject_id( ...
            data_bhv{i_cond}(i_sub).subsName)), ...
            'Subject order differs across behavior conditions at row %d.', ...
            i_sub);
    end
end
end


function id = normalize_subject_id(value)
if iscell(value),value=value{1};end
id = strtrim(char(string(value)));
end


function file = find_unique_subject_file(folder,subject_id)
items = dir(fullfile(folder,[subject_id '*']));
items = items(~[items.isdir]);
assert(isscalar(items),'Expected one file for %s in %s, found %d.', ...
    subject_id,folder,numel(items));
file = fullfile(items(1).folder,items(1).name);
end


function data = load_raw_with_rejection(file)
S = load(file);
names = fieldnames(S);
data = [];
for i = 1:numel(names)
    x = S.(names{i});
    if isstruct(x) && isfield(x,'trial') && ...
            isfield(x,'time') && isfield(x,'label')
        assert(isempty(data),'MAT file has multiple raw structures: %s',file);
        data = x;
    end
end
assert(~isempty(data),'No FieldTrip raw structure found: %s',file);
assert(isfield(data,'rejection_info') && iscell(data.rejection_info) && ...
    numel(data.rejection_info)==3, ...
    'rejection_info must contain three run cells: %s',file);
end


function [leadfields,grads] = load_three_run_leadfields(file)
S = load(file,'leadfield_individual','sensor_grad0');
assert(isfield(S,'leadfield_individual') && ...
    iscell(S.leadfield_individual) && numel(S.leadfield_individual)==3, ...
    'leadfield_individual must contain three runs: %s',file);
leadfields = S.leadfield_individual;
if isfield(S,'sensor_grad0'),grads=S.sensor_grad0;else,grads=[];end
end


function [bhv,bad_global,bad_by_run] = remove_rejected_behavior_trials( ...
    bhv,rejection_info,trials_per_run)
n_original = get_behavior_trial_count(bhv);
assert(n_original==3*trials_per_run, ...
    'Expected %d original trials, found %d.',3*trials_per_run,n_original);
bad_global = [];
bad_by_run = cell(3,1);
for i_run = 1:3
    assert(isfield(rejection_info{i_run},'bad_trials'));
    b = unique(double(rejection_info{i_run}.bad_trials(:)),'stable');
    assert(all(b>=1 & b<=trials_per_run), ...
        'Run %d rejection index outside 1:%d.',i_run,trials_per_run);
    bad_by_run{i_run} = b;
    bad_global = [bad_global;b+(i_run-1)*trials_per_run]; %#ok<AGROW>
end
bad_global = unique(bad_global,'stable');
fields = fieldnames(bhv);
for i = 1:numel(fields)
    value = bhv.(fields{i});
    if (isnumeric(value) || islogical(value) || iscell(value) || ...
            isstring(value)) && isvector(value) && numel(value)==n_original
        value(bad_global) = [];
        bhv.(fields{i}) = value;
    end
end
end


function n = get_behavior_trial_count(bhv)
candidates = {'answer1','conflict','ischange','lo_s1'};
n = [];
for i = 1:numel(candidates)
    if isfield(bhv,candidates{i})
        n = numel(bhv.(candidates{i}));
        break;
    end
end
assert(~isempty(n),'Cannot determine behavior trial count.');
end


function run_id = make_run_id_after_rejection(trials_per_run,bad_global)
run_id = repelem((1:3)',trials_per_run);
run_id(bad_global) = [];
end


function tbl = behavior_struct_to_table(bhv,idx,fields)
tbl = table();
for i = 1:numel(fields)
    f = fields{i};
    if ~isfield(bhv,f)
        if ismember(f,{'lo_chosen','lo_unchosen'}),continue;end
        error('Behavior lacks required field %s.',f);
    end
    value = bhv.(f);
    assert(numel(value)==get_behavior_trial_count(bhv), ...
        'Behavior field %s has inconsistent trial count.',f);
    tbl.(f) = double(value(idx));
    tbl.(f) = tbl.(f)(:);
end
end


function data = rmfield_if_exists(data,names)
for i = 1:numel(names)
    if isfield(data,names{i}),data=rmfield(data,names{i});end
end
end


function value = get_subject_index_vector(container,i_sub)
if iscell(container),value=container{i_sub};else,value=container(i_sub,:);end
value = double(value(:)');
value = value(isfinite(value) & value~=0);
end


function value = get_subject_scalar(container,i_sub)
if iscell(container),value=container{i_sub};else,value=container(i_sub);end
value = double(value);
assert(isscalar(value) && isfinite(value));
end


function validate_day_run_indices(idx,name)
assert(all(idx==round(idx) & idx>=1 & idx<=6), ...
    '%s must contain integers from 1 to 6.',name);
end


function [data,bhv,used_slots] = append_day_runs( ...
    all_data,all_bhv,day_slots,keep_within_day,name)
slots = day_slots(keep_within_day);
has_data = ~cellfun(@isempty,all_data(slots));
slots = slots(has_data);
assert(~isempty(slots),'%s has no usable run data.',name);
data_cells = all_data(slots);
bhv_cells = all_bhv(slots);
reference_label = data_cells{1}.label(:);
for i = 1:numel(data_cells)
    this_label = data_cells{i}.label(:);
    assert(numel(this_label)==numel(reference_label) && ...
        all(ismember(reference_label,this_label)) && ...
        all(ismember(this_label,reference_label)), ...
        ['%s retained slot %d has a different gradiometer set. ' ...
        'Check bad-channel interpolation.'],name,slots(i));
    if ~isequal(this_label,reference_label)
        cfg_reorder = [];
        cfg_reorder.channel = reference_label;
        data_cells{i} = ft_selectdata(cfg_reorder,data_cells{i});
    end
end
cfg_app = [];
cfg_app.keepsampleinfo = 'no';
data = ft_appenddata(cfg_app,data_cells{:});
bhv = vertcat(bhv_cells{:});
assert(numel(data.trial)==height(bhv), ...
    '%s appended MEG/behavior counts differ.',name);
used_slots = slots;
end


function grid = load_source_grid(file)
S = load(file);
names = fieldnames(S);
grid = [];
for i = 1:numel(names)
    x = S.(names{i});
    if isstruct(x) && isfield(x,'pos'),grid=x;break;end
end
assert(~isempty(grid),'No source grid found in %s.',file);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function result = run_sigmaV_timebeta_vector_lcmv_core( ...
    ot_data,pl_data,bhv_ot,bhv_pl,ot_lf,pl_lf,source_grid,cfg)
% RUN_SIGMAV_TIMEBETA_VECTOR_LCMV_CORE
%
% EXPLORATORY SIGMAV-BETA-DOMAIN, ORIENTATION-INVARIANT SOURCE ANALYSIS.
%
% Procedure
% ---------
% 1. Recode chosen/unchosen confidence and form two predictors:
%       SigmaV = chosen + unchosen; DeltaV = chosen - unchosen.
% 2. Fit Y(t) = b0(t) + bSigma(t)*SigmaV + bDelta(t)*DeltaV at every
%    time point on 204 signed gradients, separately for target and ITI.
% 3. Construct one SigmaV-beta covariance per day by equally averaging
%    target and ITI SigmaV-beta temporal second-order matrices.
% 4. Build one rank-2 OT/PL joint LCMV filter at every source grid point.
% 5. Split the filter into OT and PL blocks and optionally make each day
%    satisfy the same two-dimensional unit-gain constraint.
% 6. Project projection-target SigmaV beta and projection-ITI SigmaV beta
%    separately while retaining both source components.
% 7. Convert the two components to sqrt(q1^2+q2^2), average OT/PL
%    equally, and output target, ITI and target-minus-ITI source maps.
%
% Important warning
% -----------------
% The same beta time series determines the covariance and projected map
% when the construction and projection windows are the same. This is an
% effect-informed descriptive sensitivity analysis, not an independent
% confirmatory localization.
%
% This is closer to the unsigned planar-magnitude sensor analysis than a
% one-dimensional signed source orientation, but it is not algebraically
% identical to fitting a GLM after sqrt(g1^2+g2^2) sensor combination.
%
% Inputs use the same data structures as
% Only bSigma enters filter construction and source-map projection. bDelta
% is retained in sensor_beta_time and GLM diagnostics.

if nargin<8 || isempty(cfg)
    cfg = struct();
end
cfg = set_defaults(cfg);

same_beta_windows = isequal( ...
    cfg.filter_target_window,cfg.projection_target_window) && ...
    cfg.filter_iti_delay==cfg.projection_iti_delay && ...
    cfg.filter_iti_duration==cfg.projection_iti_duration;
if same_beta_windows
    warning(['filter_window and projection_window are identical. The ' ...
        'same time-resolved beta constructs C_beta and generates the ' ...
        'source maps; interpret outputs descriptively.']);
else
    warning(['Different beta windows are used for filter construction ' ...
        'and projection. The analysis remains effect-informed because ' ...
        'both windows come from the same trials/behavioral model.']);
end

validate_raw(ot_data,'ot_data');
validate_raw(pl_data,'pl_data');
bhv_ot = ensure_behavior_table(bhv_ot,numel(ot_data.trial),'OT');
bhv_pl = ensure_behavior_table(bhv_pl,numel(pl_data.trial),'PL');

%% Align each day to its representative leadfield
[ot_data,label_ot,lf_idx_ot] = align_data_leadfield( ...
    ot_data,ot_lf,cfg.channel,'OT');
[pl_data,label_pl,lf_idx_pl] = align_data_leadfield( ...
    pl_data,pl_lf,cfg.channel,'PL');
assert(numel(label_ot)==cfg.expected_n_grad, ...
    'OT has %d channels; expected %d.',numel(label_ot),cfg.expected_n_grad);
assert(numel(label_pl)==cfg.expected_n_grad, ...
    'PL has %d channels; expected %d.',numel(label_pl),cfg.expected_n_grad);

fs_ot = get_fsample(ot_data);
fs_pl = get_fsample(pl_data);
assert(abs(fs_ot-cfg.expected_fsample)<cfg.fsample_tolerance);
assert(abs(fs_pl-cfg.expected_fsample)<cfg.fsample_tolerance);

%% Independent target/ITI filter and projection cubes
[Yft_ot,time_ft_ot,info_ft_ot] = extract_fixed_window_cube( ...
    ot_data,cfg.filter_target_window,'OT filter target');
[Yft_pl,time_ft_pl,info_ft_pl] = extract_fixed_window_cube( ...
    pl_data,cfg.filter_target_window,'PL filter target');
assert_same_time(time_ft_ot,time_ft_pl,'filter target beta');

[Ypt_ot,time_pt_ot,info_pt_ot] = extract_fixed_window_cube( ...
    ot_data,cfg.projection_target_window,'OT projection target');
[Ypt_pl,time_pt_pl,info_pt_pl] = extract_fixed_window_cube( ...
    pl_data,cfg.projection_target_window,'PL projection target');
assert_same_time(time_pt_ot,time_pt_pl,'projection target beta');

[Yfi_ot,time_fi_ot,info_fi_ot] = extract_iti_window_cube( ...
    ot_data,bhv_ot,cfg.filter_iti_duration,cfg.filter_iti_delay, ...
    cfg.filter_min_iti_duration,fs_ot,'OT filter ITI');
[Yfi_pl,time_fi_pl,info_fi_pl] = extract_iti_window_cube( ...
    pl_data,bhv_pl,cfg.filter_iti_duration,cfg.filter_iti_delay, ...
    cfg.filter_min_iti_duration,fs_pl,'PL filter ITI');
assert_same_time(time_fi_ot,time_fi_pl,'filter ITI beta');

[Ypi_ot,time_pi_ot,info_pi_ot] = extract_iti_window_cube( ...
    ot_data,bhv_ot,cfg.projection_iti_duration, ...
    cfg.projection_iti_delay,cfg.projection_min_iti_duration,fs_ot, ...
    'OT projection ITI');
[Ypi_pl,time_pi_pl,info_pi_pl] = extract_iti_window_cube( ...
    pl_data,bhv_pl,cfg.projection_iti_duration, ...
    cfg.projection_iti_delay,cfg.projection_min_iti_duration,fs_pl, ...
    'PL projection ITI');
assert_same_time(time_pi_ot,time_pi_pl,'projection ITI beta');

filter_time_target = time_ft_ot(:)';
filter_time_iti = time_fi_ot(:)';
projection_time_target = time_pt_ot(:)';
projection_time_iti = time_pi_ot(:)';

%% SigmaV/DeltaV predictors and shared OT/PL scaling
[plus_ot,minus_ot,include_ot] = make_plus_minus_predictors(bhv_ot,cfg);
[plus_pl,minus_pl,include_pl] = make_plus_minus_predictors(bhv_pl,cfg);
valid_ot = include_ot & isfinite(plus_ot) & isfinite(minus_ot) & ...
    info_ft_ot.valid & info_pt_ot.valid & ...
    info_fi_ot.valid & info_pi_ot.valid & ...
    info_fi_ot.full_length & info_pi_ot.full_length;
valid_pl = include_pl & isfinite(plus_pl) & isfinite(minus_pl) & ...
    info_ft_pl.valid & info_pt_pl.valid & ...
    info_fi_pl.valid & info_pi_pl.valid & ...
    info_fi_pl.full_length & info_pi_pl.full_length;
assert(sum(valid_ot)>=cfg.min_trial, ...
    'OT has only %d valid trials.',sum(valid_ot));
assert(sum(valid_pl)>=cfg.min_trial, ...
    'PL has only %d valid trials.',sum(valid_pl));

plus_pool = [plus_ot(valid_ot);plus_pl(valid_pl)];
minus_pool = [minus_ot(valid_ot);minus_pl(valid_pl)];
scale_info = struct();
scale_info.plus_mean = mean(plus_pool,'omitnan');
scale_info.minus_mean = mean(minus_pool,'omitnan');
scale_info.plus_sd = std(plus_pool,0,'omitnan');
scale_info.minus_sd = std(minus_pool,0,'omitnan');
assert(isfinite(scale_info.plus_sd) && scale_info.plus_sd>eps, ...
    'SigmaV predictor has zero/non-finite pooled SD.');
assert(isfinite(scale_info.minus_sd) && scale_info.minus_sd>eps, ...
    'DeltaV predictor has zero/non-finite pooled SD.');
scale_info.scaling_scope = 'pooled OT/PL; each predictor scaled separately';
plus_z_ot = (plus_ot-scale_info.plus_mean)/scale_info.plus_sd;
minus_z_ot = (minus_ot-scale_info.minus_mean)/scale_info.minus_sd;
plus_z_pl = (plus_pl-scale_info.plus_mean)/scale_info.plus_sd;
minus_z_pl = (minus_pl-scale_info.minus_mean)/scale_info.minus_sd;

%% Time-resolved GLMs for all four phases
[bft_p_ot,bft_m_ot,glm_ft_ot] = fit_time_glm( ...
    Yft_ot,plus_z_ot,minus_z_ot,bhv_ot,valid_ot,cfg, ...
    'OT filter target');
[bft_p_pl,bft_m_pl,glm_ft_pl] = fit_time_glm( ...
    Yft_pl,plus_z_pl,minus_z_pl,bhv_pl,valid_pl,cfg, ...
    'PL filter target');
[bfi_p_ot,bfi_m_ot,glm_fi_ot] = fit_time_glm( ...
    Yfi_ot,plus_z_ot,minus_z_ot,bhv_ot,valid_ot,cfg, ...
    'OT filter ITI');
[bfi_p_pl,bfi_m_pl,glm_fi_pl] = fit_time_glm( ...
    Yfi_pl,plus_z_pl,minus_z_pl,bhv_pl,valid_pl,cfg, ...
    'PL filter ITI');
[bpt_p_ot,bpt_m_ot,glm_pt_ot] = fit_time_glm( ...
    Ypt_ot,plus_z_ot,minus_z_ot,bhv_ot,valid_ot,cfg, ...
    'OT projection target');
[bpt_p_pl,bpt_m_pl,glm_pt_pl] = fit_time_glm( ...
    Ypt_pl,plus_z_pl,minus_z_pl,bhv_pl,valid_pl,cfg, ...
    'PL projection target');
[bpi_p_ot,bpi_m_ot,glm_pi_ot] = fit_time_glm( ...
    Ypi_ot,plus_z_ot,minus_z_ot,bhv_ot,valid_ot,cfg, ...
    'OT projection ITI');
[bpi_p_pl,bpi_m_pl,glm_pi_pl] = fit_time_glm( ...
    Ypi_pl,plus_z_pl,minus_z_pl,bhv_pl,valid_pl,cfg, ...
    'PL projection ITI');

if cfg.beta_smooth_samples>1
    beta_cells = {bft_p_ot,bft_m_ot,bft_p_pl,bft_m_pl, ...
        bfi_p_ot,bfi_m_ot,bfi_p_pl,bfi_m_pl, ...
        bpt_p_ot,bpt_m_ot,bpt_p_pl,bpt_m_pl, ...
        bpi_p_ot,bpi_m_ot,bpi_p_pl,bpi_m_pl};
    for i_beta = 1:numel(beta_cells)
        beta_cells{i_beta} = smoothdata(beta_cells{i_beta}, ...
            2,'movmean', ...
            cfg.beta_smooth_samples,'omitnan');
    end
    [bft_p_ot,bft_m_ot,bft_p_pl,bft_m_pl, ...
        bfi_p_ot,bfi_m_ot,bfi_p_pl,bfi_m_pl, ...
        bpt_p_ot,bpt_m_ot,bpt_p_pl,bpt_m_pl, ...
        bpi_p_ot,bpi_m_ot,bpi_p_pl,bpi_m_pl] = beta_cells{:};
end

filter_target_plus_ot = bft_p_ot;
filter_target_plus_pl = bft_p_pl;
filter_iti_plus_ot = bfi_p_ot;
filter_iti_plus_pl = bfi_p_pl;
projection_target_plus_ot = bpt_p_ot;
projection_target_plus_pl = bpt_p_pl;
projection_iti_plus_ot = bpi_p_ot;
projection_iti_plus_pl = bpi_p_pl;

%% Equal target/ITI SigmaV-beta covariance per day
[C_target_ot,cov_target_ot] = beta_time_covariance( ...
    filter_target_plus_ot,cfg.beta_covariance_mode,'OT target SigmaV');
[C_iti_ot,cov_iti_ot] = beta_time_covariance( ...
    filter_iti_plus_ot,cfg.beta_covariance_mode,'OT ITI SigmaV');
[C_target_pl,cov_target_pl] = beta_time_covariance( ...
    filter_target_plus_pl,cfg.beta_covariance_mode,'PL target SigmaV');
[C_iti_pl,cov_iti_pl] = beta_time_covariance( ...
    filter_iti_plus_pl,cfg.beta_covariance_mode,'PL ITI SigmaV');
C_beta_ot = real((C_target_ot+C_iti_ot)/2);
C_beta_pl = real((C_target_pl+C_iti_pl)/2);
cov_info_ot = struct('target',cov_target_ot, ...
    'iti',cov_iti_ot,'combination','equal target/ITI mean');
cov_info_pl = struct('target',cov_target_pl, ...
    'iti',cov_iti_pl,'combination','equal target/ITI mean');
[C_ot_reg,reg_ot] = regularize_covariance( ...
    C_beta_ot,cfg.lambda,'OT beta covariance');
[C_pl_reg,reg_pl] = regularize_covariance( ...
    C_beta_pl,cfg.lambda,'PL beta covariance');
C_joint = blkdiag(C_ot_reg,C_pl_reg);
invC = pinv(C_joint);

%% Rank-2 joint source filter and separate beta projection
n_ot = numel(label_ot);
n_pl = numel(label_pl);
n_time_target = numel(projection_time_target);
n_time_iti = numel(projection_time_iti);
n_grid = numel(ot_lf.leadfield);
assert(numel(pl_lf.leadfield)==n_grid);
assert(size(source_grid.pos,1)==n_grid);
inside = get_inside_mask(source_grid,n_grid) & ...
    get_inside_mask(ot_lf,n_grid) & get_inside_mask(pl_lf,n_grid);
if isempty(cfg.grid_indices)
    grid_indices = find(inside(:))';
else
    grid_indices = intersect(find(inside(:)),cfg.grid_indices(:),'stable')';
    inside(:) = false;
    inside(grid_indices) = true;
end

source_target_time_raw_magnitude = nan(n_grid,n_time_target);
source_iti_time_raw_magnitude = nan(n_grid,n_time_iti);
source_target_time_ung_magnitude = nan(n_grid,n_time_target);
source_iti_time_ung_magnitude = nan(n_grid,n_time_iti);
source_target_time_raw_sumsq = nan(n_grid,n_time_target);
source_iti_time_raw_sumsq = nan(n_grid,n_time_iti);
source_target_time_ung_sumsq = nan(n_grid,n_time_target);
source_iti_time_ung_sumsq = nan(n_grid,n_time_iti);
gain_error_ot = nan(n_grid,1);
gain_error_pl = nan(n_grid,1);
orientation_basis_xyz = nan(n_grid,3,2,'single');
if cfg.keepfilter
    filter_ot_raw = nan(n_grid,2,n_ot,cfg.filter_storage_class);
    filter_pl_raw = nan(n_grid,2,n_pl,cfg.filter_storage_class);
    filter_ot_ung = nan(n_grid,2,n_ot,cfg.filter_storage_class);
    filter_pl_ung = nan(n_grid,2,n_pl,cfg.filter_storage_class);
else
    filter_ot_raw = [];
    filter_pl_raw = [];
    filter_ot_ung = [];
    filter_pl_ung = [];
end

for i_grid = grid_indices
    L_ot0 = ot_lf.leadfield{i_grid};
    L_pl0 = pl_lf.leadfield{i_grid};
    if isempty(L_ot0) || isempty(L_pl0)
        inside(i_grid) = false;
        continue;
    end
    L_ot = double(L_ot0(lf_idx_ot,:));
    L_pl = double(L_pl0(lf_idx_pl,:));
    L_joint = [L_ot;L_pl];

    [~,~,V] = svd(L_joint,'econ');
    r = min(cfg.source_rank,rank(L_joint));
    if r<1
        inside(i_grid) = false;
        continue;
    end
    Vr = V(:,1:r);
    if r~=2
        inside(i_grid) = false;
        continue;
    end
    Lr_ot = L_ot*Vr;
    Lr_pl = L_pl*Vr;
    Lr_joint = [Lr_ot;Lr_pl];
    G = real((Lr_joint'*invC*Lr_joint+ ...
        (Lr_joint'*invC*Lr_joint)')/2);
    if rcond(G)<cfg.min_rcond
        inside(i_grid) = false;
        continue;
    end
    W_joint = pinv(G)*(Lr_joint'*invC); % 2 x 408
    W_ot = real(W_joint(:,1:n_ot));
    W_pl = real(W_joint(:,n_ot+1:end));

    % Make each day separately satisfy W_day*L_day=I2.
    G_ot = W_ot*Lr_ot;
    G_pl = W_pl*Lr_pl;
    if cfg.apply_day_gain_calibration
        if rcond(G_ot)<cfg.min_rcond || rcond(G_pl)<cfg.min_rcond
            inside(i_grid) = false;
            continue;
        end
        W_ot = pinv(G_ot)*W_ot;
        W_pl = pinv(G_pl)*W_pl;
    end
    gain_error_ot(i_grid) = norm(W_ot*Lr_ot-eye(2),'fro');
    gain_error_pl(i_grid) = norm(W_pl*Lr_pl-eye(2),'fro');

    % Two-component unit white-noise-gain normalization. There is no
    % independent sensor-noise covariance in this beta-domain analysis.
    [A_ot,ok_ot] = invsqrt_psd_local(W_ot*W_ot',cfg.min_rcond);
    [A_pl,ok_pl] = invsqrt_psd_local(W_pl*W_pl',cfg.min_rcond);
    if ~ok_ot || ~ok_pl
        inside(i_grid) = false;
        continue;
    end
    W_ot_ung = A_ot*W_ot;
    W_pl_ung = A_pl*W_pl;

    [source_target_time_raw_magnitude(i_grid,:), ...
        source_target_time_raw_sumsq(i_grid,:)] = ...
        day_balanced_vector_measures( ...
        W_ot,projection_target_plus_ot,W_pl,projection_target_plus_pl);
    [source_iti_time_raw_magnitude(i_grid,:), ...
        source_iti_time_raw_sumsq(i_grid,:)] = ...
        day_balanced_vector_measures( ...
        W_ot,projection_iti_plus_ot,W_pl,projection_iti_plus_pl);
    [source_target_time_ung_magnitude(i_grid,:), ...
        source_target_time_ung_sumsq(i_grid,:)] = ...
        day_balanced_vector_measures( ...
        W_ot_ung,projection_target_plus_ot, ...
        W_pl_ung,projection_target_plus_pl);
    [source_iti_time_ung_magnitude(i_grid,:), ...
        source_iti_time_ung_sumsq(i_grid,:)] = ...
        day_balanced_vector_measures( ...
        W_ot_ung,projection_iti_plus_ot, ...
        W_pl_ung,projection_iti_plus_pl);
    orientation_basis_xyz(i_grid,:,:) = reshape(single(Vr),1,3,2);
    if cfg.keepfilter
        filter_ot_raw(i_grid,:,:) = reshape( ...
            cast(W_ot,cfg.filter_storage_class),1,2,n_ot);
        filter_pl_raw(i_grid,:,:) = reshape( ...
            cast(W_pl,cfg.filter_storage_class),1,2,n_pl);
        filter_ot_ung(i_grid,:,:) = reshape( ...
            cast(W_ot_ung,cfg.filter_storage_class),1,2,n_ot);
        filter_pl_ung(i_grid,:,:) = reshape( ...
            cast(W_pl_ung,cfg.filter_storage_class),1,2,n_pl);
    end
end

valid_grid = inside & ...
    all(isfinite(source_target_time_raw_magnitude),2) & ...
    all(isfinite(source_iti_time_raw_magnitude),2) & ...
    all(isfinite(source_target_time_ung_magnitude),2) & ...
    all(isfinite(source_iti_time_ung_magnitude),2);
assert(any(valid_grid),'No valid source filters were produced.');

maps_magnitude = summarize_vector_measure( ...
    source_grid,valid_grid, ...
    source_target_time_raw_magnitude,source_iti_time_raw_magnitude, ...
    source_target_time_ung_magnitude,source_iti_time_ung_magnitude, ...
    'magnitude');
maps_sumsq = summarize_vector_measure( ...
    source_grid,valid_grid, ...
    source_target_time_raw_sumsq,source_iti_time_raw_sumsq, ...
    source_target_time_ung_sumsq,source_iti_time_ung_sumsq, ...
    'sumsq');

%% Output
result = struct();
result.method = [ ...
    'time-resolved target/ITI SigmaV beta -> equal target/ITI C_beta/day ' ...
    '-> rank-2 blockdiag joint LCMV -> separate target/ITI SigmaV-beta ' ...
    'projection -> two-component vector magnitude and/or sum of squares'];
result.analysis_warning = [ ...
    'Effect-informed descriptive analysis; target and ITI beta from the ' ...
    'same trials/behavioral model construct the filter and source maps.'];
result.maps_magnitude = maps_magnitude;
result.maps_sumsq = maps_sumsq;
switch cfg.vector_output_mode
    case 'magnitude'
        result.maps = maps_magnitude;
        result.dynamic = pack_dynamic_measure( ...
            source_target_time_raw_magnitude,source_iti_time_raw_magnitude, ...
            source_target_time_ung_magnitude,source_iti_time_ung_magnitude);
    case 'sumsq'
        result.maps = maps_sumsq;
        result.dynamic = pack_dynamic_measure( ...
            source_target_time_raw_sumsq,source_iti_time_raw_sumsq, ...
            source_target_time_ung_sumsq,source_iti_time_ung_sumsq);
    case 'both'
        % Backward-compatible aliases remain magnitude in 'both' mode.
        result.maps = maps_magnitude;
        result.dynamic = pack_dynamic_measure( ...
            source_target_time_raw_magnitude,source_iti_time_raw_magnitude, ...
            source_target_time_ung_magnitude,source_iti_time_ung_magnitude);
        result.dynamic_magnitude = result.dynamic;
        result.dynamic_sumsq = pack_dynamic_measure( ...
            source_target_time_raw_sumsq,source_iti_time_raw_sumsq, ...
            source_target_time_ung_sumsq,source_iti_time_ung_sumsq);
    otherwise
        error('Unknown vector_output_mode: %s',cfg.vector_output_mode);
end
result.time = struct('filter_target',filter_time_target, ...
    'filter_iti',filter_time_iti, ...
    'projection_target',projection_time_target, ...
    'projection_iti',projection_time_iti);
result.sensor_beta_time = struct();
result.sensor_beta_time.filter_target = pack_plus_minus_beta( ...
    bft_p_ot,bft_m_ot,bft_p_pl,bft_m_pl);
result.sensor_beta_time.filter_iti = pack_plus_minus_beta( ...
    bfi_p_ot,bfi_m_ot,bfi_p_pl,bfi_m_pl);
result.sensor_beta_time.projection_target = pack_plus_minus_beta( ...
    bpt_p_ot,bpt_m_ot,bpt_p_pl,bpt_m_pl);
result.sensor_beta_time.projection_iti = pack_plus_minus_beta( ...
    bpi_p_ot,bpi_m_ot,bpi_p_pl,bpi_m_pl);
result.predictor_scale = scale_info;
result.glm = struct('filter_target_ot',glm_ft_ot, ...
    'filter_target_pl',glm_ft_pl,'filter_iti_ot',glm_fi_ot, ...
    'filter_iti_pl',glm_fi_pl,'projection_target_ot',glm_pt_ot, ...
    'projection_target_pl',glm_pt_pl, ...
    'projection_iti_ot',glm_pi_ot,'projection_iti_pl',glm_pi_pl);
result.valid_trial_ot = valid_ot;
result.valid_trial_pl = valid_pl;
result.window_info = struct('filter_target_ot',info_ft_ot, ...
    'filter_target_pl',info_ft_pl,'filter_iti_ot',info_fi_ot, ...
    'filter_iti_pl',info_fi_pl,'projection_target_ot',info_pt_ot, ...
    'projection_target_pl',info_pt_pl, ...
    'projection_iti_ot',info_pi_ot,'projection_iti_pl',info_pi_pl);
result.label_ot = label_ot;
result.label_pl = label_pl;
result.beta_covariance = struct( ...
    'C_ot',C_beta_ot,'C_pl',C_beta_pl, ...
    'C_target_ot',C_target_ot,'C_iti_ot',C_iti_ot, ...
    'C_target_pl',C_target_pl,'C_iti_pl',C_iti_pl, ...
    'C_ot_regularized',C_ot_reg,'C_pl_regularized',C_pl_reg, ...
    'C_joint_regularized',C_joint, ...
    'info_ot',cov_info_ot,'info_pl',cov_info_pl, ...
    'regularization_ot',reg_ot,'regularization_pl',reg_pl);
result.orientation_basis_xyz = orientation_basis_xyz;
result.filter = struct('ot_raw',filter_ot_raw,'pl_raw',filter_pl_raw, ...
    'ot_unitnoisegain',filter_ot_ung, ...
    'pl_unitnoisegain',filter_pl_ung);
result.quality = struct( ...
    'n_valid_grid',sum(valid_grid), ...
    'max_unit_gain_error_ot',max(gain_error_ot(valid_grid), ...
        [],'omitnan'), ...
    'max_unit_gain_error_pl',max(gain_error_pl(valid_grid), ...
        [],'omitnan'), ...
    'rank_C_beta_ot',rank(C_beta_ot), ...
    'rank_C_beta_pl',rank(C_beta_pl), ...
    'rcond_C_ot_regularized',rcond(C_ot_reg), ...
    'rcond_C_pl_regularized',rcond(C_pl_reg));
result.cfg = cfg;
end


function cfg = set_defaults(cfg)
if isfield(cfg,'target_window') && ~isempty(cfg.target_window)
    legacy_target_window = cfg.target_window;
else
    legacy_target_window = [-0.5 -0.3];
end
defaults = struct();
defaults.filter_target_window = legacy_target_window;
defaults.projection_target_window = legacy_target_window;
defaults.filter_iti_delay = 0.3;
defaults.filter_iti_duration = 0.2;
defaults.filter_min_iti_duration = 0.2;
defaults.projection_iti_delay = 0.3;
defaults.projection_iti_duration = 0.2;
defaults.projection_min_iti_duration = 0.2;
defaults.condition_include = 1:4;
defaults.require_conflict = true;
defaults.channel = 'MEGGRAD';
defaults.expected_n_grad = 204;
defaults.expected_fsample = 250;
defaults.fsample_tolerance = 0.5;
defaults.lambda = 0.10;
% secondmoment retains a sustained beta effect; covariance removes each
% channel's temporal mean and is recommended only as a sensitivity check.
defaults.beta_covariance_mode = 'secondmoment'; % covariance | secondmoment
defaults.beta_smooth_samples = 1;
defaults.source_rank = 2;
defaults.min_trial = 20;
% Keep these false to match the current sensor-level chosen/unchosen GLM.
defaults.add_change_nuisance = false;
defaults.add_condition_nuisance = false;
defaults.add_run_nuisance = false;
defaults.apply_day_gain_calibration = true;
defaults.vector_output_mode = 'magnitude';
defaults.min_rcond = 1e-10;
defaults.keepfilter = false;
defaults.filter_storage_class = 'single';
defaults.grid_indices = [];
names = fieldnames(defaults);
for i = 1:numel(names)
    if ~isfield(cfg,names{i}) || isempty(cfg.(names{i}))
        cfg.(names{i}) = defaults.(names{i});
    end
end
cfg.beta_covariance_mode = char(lower(string( ...
    cfg.beta_covariance_mode)));
cfg.vector_output_mode = char(lower(string(cfg.vector_output_mode)));
assert(ismember(cfg.beta_covariance_mode, ...
    {'covariance','secondmoment'}));
assert(ismember(cfg.vector_output_mode, ...
    {'magnitude','sumsq','both'}), ...
    'vector_output_mode must be magnitude, sumsq, or both.');
assert(cfg.lambda>0,'Positive ridge regularization is required.');
assert(cfg.beta_smooth_samples>=1 && ...
    mod(cfg.beta_smooth_samples,1)==0);
assert(cfg.source_rank==2, ...
    'This vector-magnitude implementation requires source_rank=2.');
assert(cfg.min_rcond>0);
assert(ismember(cfg.filter_storage_class,{'single','double'}));
assert(numel(cfg.filter_target_window)==2 && ...
    cfg.filter_target_window(1)<cfg.filter_target_window(2));
assert(numel(cfg.projection_target_window)==2 && ...
    cfg.projection_target_window(1)<cfg.projection_target_window(2));
assert(cfg.filter_iti_duration>=cfg.filter_min_iti_duration && ...
    cfg.filter_min_iti_duration>0);
assert(cfg.projection_iti_duration>=cfg.projection_min_iti_duration && ...
    cfg.projection_min_iti_duration>0);
end


function validate_raw(data,name)
assert(isstruct(data) && isfield(data,'trial') && ...
    isfield(data,'time') && isfield(data,'label'), ...
    '%s is not a FieldTrip raw structure.',name);
assert(numel(data.trial)==numel(data.time));
end


function bhv = ensure_behavior_table(bhv,n_trial,name)
required = {'conflict','ischange','lo_s1','lo_o1','condition','run_id', ...
    'oiti','osti2','rtime2'};
if ~istable(bhv)
    assert(isstruct(bhv) && isscalar(bhv));
    tbl = table();
    for i = 1:numel(required)
        f = required{i};
        assert(isfield(bhv,f),'%s behavior lacks %s.',name,f);
        tbl.(f) = double(bhv.(f)(:));
    end
    if isfield(bhv,'run_uid')
        tbl.run_uid = double(bhv.run_uid(:));
    end
    bhv = tbl;
end
for i = 1:numel(required)
    assert(ismember(required{i},bhv.Properties.VariableNames), ...
        '%s behavior lacks %s.',name,required{i});
end
assert(height(bhv)==n_trial, ...
    '%s behavior/MEG trial counts differ.',name);
if ~ismember('run_uid',bhv.Properties.VariableNames)
    bhv.run_uid = (double(bhv.condition)-1)*3+double(bhv.run_id);
end
end


function [data,labels,lf_idx] = align_data_leadfield( ...
    data,leadfield,channel_selection,name)
cfg_sel = [];
cfg_sel.channel = channel_selection;
data = ft_selectdata(cfg_sel,data);
keep = ismember(data.label,leadfield.label);
labels = data.label(keep);
assert(~isempty(labels),'%s has no channels shared with leadfield.',name);
cfg_sel = [];
cfg_sel.channel = labels;
data = ft_selectdata(cfg_sel,data);
[found,lf_idx] = ismember(data.label,leadfield.label);
assert(all(found));
labels = data.label;
end


function fs = get_fsample(data)
if isfield(data,'fsample') && isfinite(data.fsample)
    fs = double(data.fsample);
else
    fs = 1/median(diff(double(data.time{1})));
end
end


function [Y,time_ref,info] = extract_fixed_window_cube(data,window,name)
n_trial = numel(data.trial);
n_chan = numel(data.label);
t0 = double(data.time{2}(:)');
idx0 = find(t0>=window(1) & t0<window(2));
assert(numel(idx0)>=2,'%s contains fewer than two samples.',name);
time_ref = t0(idx0);
n_time = numel(time_ref);
% fs = get_fsample(data);
% sample_start = round(window(1)*fs);
% sample_end   = round(window(2)*fs)-1;
% time_ref = (sample_start:sample_end) ./ fs;
% n_time = numel(time_ref);

Y = nan(n_trial,n_chan,n_time);
valid = false(n_trial,1);

for i_trial = 1:n_trial
    t = double(data.time{i_trial}(:)');
    idx = find(t>=window(1) & t<window(2));
    if numel(idx)~=n_time || max(abs(t(idx)-time_ref))>1e-8
        continue;
    end
    x = double(data.trial{i_trial}(:,idx));
    if ~all(isfinite(x(:)))
        continue;
    end
    Y(i_trial,:,:) = x;
    valid(i_trial) = true;
end

info = struct('name',name,'window',window,'time',time_ref, ...
    'n_time',n_time,'valid',valid,'n_valid',sum(valid), ...
    'n_invalid',sum(~valid));
end


function [Y,relative_time,info] = extract_iti_window_cube( ...
    data,bhv,duration,delay,min_duration,fs,name)
n_trial = numel(data.trial);
n_chan = numel(data.label);
n_requested = round(duration*fs);
n_min = ceil(min_duration*fs-1e-8);
assert(n_requested>=2 && n_min>=2 && n_min<=n_requested);
relative_time = (0:n_requested-1)/fs;
Y = nan(n_trial,n_chan,n_requested);
valid = false(n_trial,1);
full_length = false(n_trial,1);
n_available = zeros(n_trial,1);
start_sec = double(bhv.oiti)-double(bhv.osti2)- ...
    double(bhv.rtime2)+delay;

for i_trial = 1:n_trial
    t = double(data.time{i_trial}(:)');
    if ~isfinite(start_sec(i_trial)) || isempty(t),continue;end
    % idx_start = round((start_sec(i_trial)-t(1))*fs)+1;
    idx_start = round((start_sec(i_trial)-t(1))*fs)+1;
    if idx_start<1 || idx_start>numel(t),continue;end
    time_error = abs(t(idx_start)-start_sec(i_trial));
    if time_error > 0.51/fs;continue;end   
    n_take = min(n_requested,numel(t)-idx_start+1);
    if n_take<n_min,continue;end
    idx = idx_start+(0:n_take-1);
    X = double(data.trial{i_trial}(:,idx));
    if ~all(isfinite(X(:))),continue;end
    Y(i_trial,:,1:n_take) = X;
    valid(i_trial) = true;
    full_length(i_trial) = n_take==n_requested;
    n_available(i_trial) = n_take;
end
info = struct('name',name,'duration',duration,'delay',delay, ...
    'min_duration',min_duration,'time',relative_time, ...
    'start_sec',start_sec,'valid',valid,'full_length',full_length, ...
    'n_valid',sum(valid),'n_full_length',sum(full_length), ...
    'n_available_sample',n_available);
end


function assert_same_time(a,b,name)
a = double(a(:)');
b = double(b(:)');
assert(numel(a)==numel(b) && max(abs(a-b))<1e-8, ...
    '%s OT/PL time axes differ.',name);
end


function [plus,minus,include] = make_plus_minus_predictors(bhv,cfg)
change = double(bhv.ischange)==1;
self_conf = double(bhv.lo_s1);
other_conf = double(bhv.lo_o1);
chosen = self_conf;
unchosen = other_conf;
chosen(change) = other_conf(change);
unchosen(change) = self_conf(change);

% Guard against an inconsistent behavioral recoding. The supplied
% lo_chosen/lo_unchosen columns must agree with self/other + ischange.
if all(ismember({'lo_chosen','lo_unchosen'},bhv.Properties.VariableNames))
    chosen_saved = double(bhv.lo_chosen);
    unchosen_saved = double(bhv.lo_unchosen);
    check_idx = isfinite(chosen) & isfinite(unchosen) & ...
        isfinite(chosen_saved) & isfinite(unchosen_saved);
    assert(all(abs(chosen(check_idx)-chosen_saved(check_idx))<1e-10) && ...
        all(abs(unchosen(check_idx)-unchosen_saved(check_idx))<1e-10), ...
        'lo_chosen/lo_unchosen disagree with self/other recoding.');
end

plus = chosen + unchosen;   % SigmaV: total confidence/evidence level
minus = chosen - unchosen; % DeltaV: chosen bias/evidence difference

include = ismember(double(bhv.condition),cfg.condition_include);
if cfg.require_conflict
    include = include & double(bhv.conflict)==1;
end
end


function [plus_beta,minus_beta,info] = fit_time_glm( ...
    Y,plus_z,minus_z,bhv,valid,cfg,name)
valid = logical(valid(:));
n = sum(valid);
D = [ones(n,1),plus_z(valid),minus_z(valid)];
column_names = {'intercept','SigmaV_plus','DeltaV_minus'};
assert(rank(D)==size(D,2), ...
    '%s primary design is rank deficient.',name);

candidate = zeros(n,0);
candidate_names = {};
if cfg.add_change_nuisance
    candidate = [candidate,double(bhv.ischange(valid))];
    candidate_names{end+1} = 'ischange';
end
if cfg.add_condition_nuisance
    [D0,N0] = make_dummies(double(bhv.condition(valid)),'condition');
    candidate = [candidate,D0];
    candidate_names = [candidate_names,N0];
end
if cfg.add_run_nuisance
    [D0,N0] = make_dummies(double(bhv.run_uid(valid)),'run_uid');
    candidate = [candidate,D0];
    candidate_names = [candidate_names,N0];
end

dropped = {};
for i = 1:size(candidate,2)
    if rank([D,candidate(:,i)])>rank(D)
        D = [D,candidate(:,i)]; %#ok<AGROW>
        column_names{end+1} = candidate_names{i}; %#ok<AGROW>
    else
        dropped{end+1} = candidate_names{i}; %#ok<AGROW>
    end
end

n_chan = size(Y,2);
n_time = size(Y,3);
Y2 = reshape(Y(valid,:,:),n,n_chan*n_time);
B = D\Y2;
plus_beta = reshape(B(2,:),n_chan,n_time);
minus_beta = reshape(B(3,:),n_chan,n_time);
info = struct('day',name,'n_trial',n,'design_columns',{column_names}, ...
    'dropped_redundant_nuisance',{dropped},'rank',rank(D), ...
    'n_column',size(D,2),'Y_standardized',false, ...
    'primary_model','Y ~ SigmaV_plus + DeltaV_minus', ...
    'plus_minus_correlation',corr(plus_z(valid),minus_z(valid)));
end


function [D,names] = make_dummies(code,prefix)
levels = unique(code(:),'stable');
D = zeros(numel(code),max(numel(levels)-1,0));
names = cell(1,size(D,2));
for i = 2:numel(levels)
    D(:,i-1) = code==levels(i);
    names{i-1} = sprintf('%s_%g',prefix,levels(i));
end
end


function [C,info] = beta_time_covariance(B,mode,name)
assert(size(B,2)>=2,'%s beta has fewer than two time samples.',name);
assert(all(isfinite(B(:))),'%s beta contains non-finite values.',name);
switch mode
    case 'covariance'
        B_used = B-mean(B,2);
        denominator = size(B,2)-1;
    case 'secondmoment'
        B_used = B;
        denominator = size(B,2);
    otherwise
        error('Unknown beta covariance mode: %s',mode);
end
C = (B_used*B_used')/denominator;
C = real((C+C')/2);
info = struct('name',name,'mode',mode,'n_time',size(B,2), ...
    'rank',rank(C),'rcond',rcond(C),'trace',trace(C));
end


function [magnitude,sumsq] = day_balanced_vector_measures( ...
    W_ot,B_ot,W_pl,B_pl)
Q_ot = W_ot*B_ot;
Q_pl = W_pl*B_pl;
assert(size(Q_ot,1)==2 && size(Q_pl,1)==2 && ...
    size(Q_ot,2)==size(Q_pl,2));
S_ot = sum(abs(Q_ot).^2,1);
S_pl = sum(abs(Q_pl).^2,1);
sumsq = (S_ot+S_pl)/2;
magnitude = (sqrt(S_ot)+sqrt(S_pl))/2;
end


function maps = summarize_vector_measure( ...
    source_grid,valid_grid,target_raw_time,iti_raw_time, ...
    target_ung_time,iti_ung_time,measure_name)
target_raw = mean(target_raw_time,2,'omitnan');
iti_raw = mean(iti_raw_time,2,'omitnan');
target_ung = mean(target_ung_time,2,'omitnan');
iti_ung = mean(iti_ung_time,2,'omitnan');
target_minus_iti_raw = target_raw-iti_raw;
target_minus_iti_ung = target_ung-iti_ung;
target_raw(~valid_grid) = NaN;
iti_raw(~valid_grid) = NaN;
target_ung(~valid_grid) = NaN;
iti_ung(~valid_grid) = NaN;
target_minus_iti_raw(~valid_grid) = NaN;
target_minus_iti_ung(~valid_grid) = NaN;
maps = struct();
maps.target_raw = make_effect_source(source_grid,valid_grid,target_raw, ...
    ['target_sigmaV_' measure_name]);
maps.iti_raw = make_effect_source(source_grid,valid_grid,iti_raw, ...
    ['iti_sigmaV_' measure_name]);
maps.target_minus_iti_raw = make_effect_source( ...
    source_grid,valid_grid,target_minus_iti_raw, ...
    ['target_minus_iti_sigmaV_' measure_name]);
maps.target_unitnoisegain = make_effect_source( ...
    source_grid,valid_grid,target_ung,['target_sigmaV_' measure_name]);
maps.iti_unitnoisegain = make_effect_source( ...
    source_grid,valid_grid,iti_ung,['iti_sigmaV_' measure_name]);
maps.target_minus_iti_unitnoisegain = make_effect_source( ...
    source_grid,valid_grid,target_minus_iti_ung, ...
    ['target_minus_iti_sigmaV_' measure_name]);
end


function dynamic = pack_dynamic_measure(target_raw,iti_raw,target_ung,iti_ung)
dynamic = struct('target_raw',target_raw,'iti_raw',iti_raw, ...
    'target_unitnoisegain',target_ung,'iti_unitnoisegain',iti_ung);
end


function out = pack_plus_minus_beta(plus_ot,minus_ot,plus_pl,minus_pl)
out = struct('plus_ot',plus_ot,'minus_ot',minus_ot, ...
    'plus_pl',plus_pl,'minus_pl',minus_pl, ...
    'plus_definition','chosen + unchosen', ...
    'minus_definition','chosen - unchosen');
end


function [A,ok] = invsqrt_psd_local(M,min_rcond)
M = real((M+M')/2);
[V,D] = eig(M,'vector');
D = real(D);
dmax = max(D);
ok = all(isfinite(D)) && dmax>0 && min(D)>dmax*min_rcond;
if ok
    A = V*diag(1./sqrt(D))*V';
    A = real((A+A')/2);
else
    A = nan(size(M));
end
end


function [Creg,info] = regularize_covariance(C,lambda,name)
n = size(C,1);
scale = trace(C)/n;
assert(isfinite(scale) && scale>0, ...
    '%s has zero or invalid variance.',name);
Creg = C+lambda*scale*eye(n);
Creg = real((Creg+Creg')/2);
info = struct('name',name,'lambda',lambda,'scale',scale, ...
    'rank_before',rank(C),'rank_after',rank(Creg), ...
    'rcond_before',rcond(C),'rcond_after',rcond(Creg));
end


function source = make_effect_source(grid,inside,effect,effect_name)
source = struct();
source.pos = grid.pos;
source.inside = logical(inside(:));
source.dimord = 'pos';
if isfield(grid,'dim'),source.dim = grid.dim;end
if isfield(grid,'unit'),source.unit = grid.unit;end
if isfield(grid,'coordsys'),source.coordsys = grid.coordsys;end
source.effect = double(effect(:));
source.effect_name = effect_name;
source.pow = source.effect;
end


function inside = get_inside_mask(grid,n_grid)
inside = true(n_grid,1);
if ~isfield(grid,'inside'),return;end
if islogical(grid.inside)
    inside = grid.inside(:);
else
    inside = false(n_grid,1);
    inside(double(grid.inside(:))) = true;
end
assert(numel(inside)==n_grid);
end
