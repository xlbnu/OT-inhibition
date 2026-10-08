%% basic config set
cfg_ppm = [];
cfg_ppm.time_window = [-3 1.5];
cfg_ppm.fsample = 250;

% Do not duplicate frequency points between low and high bands, Avoid ft_appendfreq Two windows appear after 36 Hz.
cfg_ppm.low_foi = 2:1:30;
cfg_ppm.high_foi = 30:2:120;
cfg_ppm.toi_step = 0.05;

cfg_ppm.low_t_ftimwin = 0.4;
cfg_ppm.high_t_ftimwin = 0.25;
cfg_ppm.high_tapsmofrq = 10;

cfg_ppm.trim_percent = 4;
cfg_ppm.trials_per_run = 80;

base_path2 = fullfile( ...
    'F:\xianliang\exp_data\OT_N2_MEG\sensor_level_TFR', ...
    sprintf('pressPutton102(%dhz)', cfg_ppm.fsample));

input_root = 'F:\xianliang\exp_data\OT_N2_MEG\MEG-Preproocess2\pressPutton';
base_path = {
    fullfile(input_root, sprintf('social_OT_306(%dhz)', cfg_ppm.fsample))
    fullfile(input_root, sprintf('social_PL_306(%dhz)', cfg_ppm.fsample))
    fullfile(input_root, sprintf('nonsocial_OT_306(%dhz)', cfg_ppm.fsample))
    fullfile(input_root, sprintf('nonsocial_PL_306(%dhz)', cfg_ppm.fsample))
    };

load('E:\xianliang\matlab_m\social_decision_m\data\OT_all\BHV_MEG_MRS_data_all\bhv_raw_data\bhv_data_MEG_19.mat');
base_data = {ots_data; pls_data; otn_data; pln_data};
cfg_ppm.trim_percent = 4;

%% conflict-congruent time-frequency
[tfr_merged,tfr_by_condition] = response_locked_tfr(base_path, base_path2, base_data, 'conflict_congruent', cfg_ppm);


%% conflict-congruent time-frequency
[tfr_merged,tfr_by_condition] = response_locked_tfr(base_path, base_path2, base_data, 'left_right', cfg_ppm);












function [tfr_merged,tfr_by_condition] = response_locked_tfr( ...
    base_path,base_path2,base_data,grouping,cfg_ppm)
% RESPONSE_LOCKED_TFR Response-locked MEG power for a selected trial grouping.
% Inputs follow the original analysis flow:
% base_path: four input directories; [] uses the original default paths.
% base_path2: output directory; [] uses response_tfr_results beside this file.
% base_data: {ots_meg;pls_meg;otn_meg;pln_meg}, in this fixed order.
% grouping: 'conflict_congruent' or 'left_right' (required).
% cfg_ppm: optional spectral/trim settings; defaults are supplied below.
% No realignment is performed: MEG input time zero must already be response.
% Each call returns and saves two structures: merged and per-condition results.
narginchk(4,5);
if nargin<5 || isempty(cfg_ppm), cfg_ppm=struct; end
assert(isstruct(cfg_ppm) && isscalar(cfg_ppm),'cfg_ppm must be a scalar struct.');
assert(ischar(grouping) || (isstring(grouping) && isscalar(grouping)), ...
    'response_tfr:InvalidGrouping','grouping must be a single text value.');
grouping=char(grouping);
assert(ismember(grouping,{'conflict_congruent','left_right'}), ...
    'response_tfr:InvalidGrouping','Use conflict_congruent or left_right.');
cfg_ppm.grouping=grouping;
cfg_ppm = configure_analysis(cfg_ppm);
if isempty(base_path)
    input_root='F:\xianliang\exp_data\OT_N2_MEG\MEG-Preproocess2\pressPutton';
    source_dirs={'social_OT_306','social_PL_306','nonsocial_OT_306','nonsocial_PL_306'};
    base_path=cellfun(@(s) fullfile(input_root,sprintf('%s(%dhz)',s,cfg_ppm.fsample)), ...
        source_dirs,'UniformOutput',false);
end
if isempty(base_path2)
    base_path2=fullfile(fileparts(mfilename('fullpath')),'response_tfr_results');
end
[tfr_merged,tfr_by_condition] = compute_response_power(base_path,base_data,cfg_ppm);
if ~isfolder(base_path2), mkdir(base_path2); end
merged_file=fullfile(base_path2,sprintf('TFR_merged_%s_%dhz.mat',cfg_ppm.grouping,cfg_ppm.fsample));
condition_file=fullfile(base_path2,sprintf('TFR_by_condition_%s_%dhz.mat',cfg_ppm.grouping,cfg_ppm.fsample));
save(merged_file,'tfr_merged','-v7.3');
save(condition_file,'tfr_by_condition','-v7.3');
fprintf('Saved merged results: %s\nSaved condition results: %s\n',merged_file,condition_file);
end

function cfg = configure_analysis(cfg)
defaults=struct('time_window',[-3 1.5],'fsample',250, ...
    'low_foi',2:35,'high_foi',36:2:120,'toi_step',0.05, ...
    'low_t_ftimwin',0.4,'high_t_ftimwin',0.25,'high_tapsmofrq',6, ...
    'trim_percent',4,'trials_per_run',80,'min_trials',1, ...
    'condition_names',{{'social_OT','social_PL','nonsocial_OT','nonsocial_PL'}});
for f=fieldnames(defaults)'
    if ~isfield(cfg,f{1}), cfg.(f{1})=defaults.(f{1}); end
end
cfg.grouping=char(cfg.grouping);
assert(ismember(cfg.grouping,{'conflict_congruent','left_right'}),'Invalid grouping.');
assert(isnumeric(cfg.time_window) && numel(cfg.time_window)==2 && ...
    all(isfinite(cfg.time_window)) && diff(cfg.time_window)>0,'Invalid time window.');
assert(cfg.fsample>0 && isfinite(cfg.fsample) && cfg.toi_step>0,'Invalid sampling/time step.');
assert(isnumeric(cfg.trim_percent) && isscalar(cfg.trim_percent) && ...
    isfinite(cfg.trim_percent) && cfg.trim_percent>=0 && cfg.trim_percent<100, ...
    'trim_percent must be in [0,100); 0 disables trimming.');
assert(cfg.min_trials>=1 && cfg.min_trials==fix(cfg.min_trials),'Invalid min_trials.');
assert(cfg.trials_per_run>=1 && cfg.trials_per_run==fix(cfg.trials_per_run),'Invalid run size.');
assert(all(diff(cfg.low_foi)>0) && all(diff(cfg.high_foi)>0) && ...
    cfg.low_foi(1)==2 && cfg.high_foi(end)==120 && ...
    max(cfg.low_foi)<min(cfg.high_foi),'Frequency bands must be ordered and disjoint.');
assert(max([cfg.low_foi,cfg.high_foi])<cfg.fsample/2,'Frequencies exceed Nyquist limit.');
end

function [merged,byCondition] = compute_response_power(paths,behavior,cfg)
assert(numel(paths)==4 && numel(behavior)==4 && numel(cfg.condition_names)==4, ...
    'Expected four experimental conditions.');
assert(numel(unique(cfg.condition_names))==4 && ...
    all(cellfun(@isvarname,cfg.condition_names)),'Invalid/duplicate condition names.');
n=numel(behavior{1}); assert(n>0,'No subjects.');
subjectID=cell(n,1);
for g=1:4
    assert(numel(behavior{g})==n,'Subject counts differ.');
    for i=1:n
        id=behavior{g}(i).subsName;
        if iscell(id), id=id{1}; end
        id=char(string(id));
        if g==1, subjectID{i}=id; else, assert(strcmp(id,subjectID{i}),'Subject order differs.'); end
    end
end
assert(numel(unique(subjectID))==n,'Duplicate subject IDs.');
switch cfg.grouping
    case 'conflict_congruent', groups={'conflict','congruent'};
    case 'left_right', groups={'left','right'};
end
merged=struct('settings',cfg,'group_names',{groups},'all_subjects',{subjectID});
byCondition=merged;
for k=1:numel(groups)
    merged.(groups{k})=struct('subjects',{cell(n,1)},'trial_log',{cell(n,1)});
end
for g=1:4
    for k=1:numel(groups)
        byCondition.(cfg.condition_names{g}).(groups{k})= ...
            struct('subjects',{cell(n,1)},'trial_log',{cell(n,1)});
    end
end
for i=1:n
    fprintf('Response power: subject %d/%d (%s)\n',i,n,subjectID{i});
    spectra=cell(1,4); bhv=cell(1,4);
    for g=1:4
        files=dir(fullfile(paths{g},[subjectID{i},'*.mat']));
        assert(isscalar(files),'Subject %s condition %d: expected exactly one input file.',subjectID{i},g);
        loaded=load(fullfile(files.folder,files.name),'data_combined_pre');
        assert(isfield(loaded,'data_combined_pre'),'Missing data_combined_pre.');
        raw=loaded.data_combined_pre;
        assert(abs(raw.fsample-cfg.fsample)<1e-6,'Input fsample does not match configuration.');
        [bhv{g},raw]=align_trials(behavior{g}(i),raw,cfg,g);
        spectra{g}=compute_spectrum(raw,cfg);
        assert(size(spectra{g}.powspctrm,1)==numel(bhv{g}.answer2),'Trial alignment changed.');
        for k=1:numel(groups)
            [average,log]=select_and_average(spectra{g},bhv{g},groups{k},cfg);
            byCondition.(cfg.condition_names{g}).(groups{k}).subjects{i}=average;
            byCondition.(cfg.condition_names{g}).(groups{k}).trial_log{i}=log;
        end
        clear raw loaded
    end
    % Pool trial powers, not condition averages: trials have equal weight.
    acfg=struct('parameter','powspctrm','appenddim','rpt');
    pooled=ft_appendfreq(acfg,spectra{:});
    combined=struct;
    for f=fieldnames(bhv{1})'
        pieces=cellfun(@(b) b.(f{1}),bhv,'UniformOutput',false);
        combined.(f{1})=vertcat(pieces{:});
    end
    assert(size(pooled.powspctrm,1)==numel(combined.answer2),'Pooled trial order mismatch.');
    for k=1:numel(groups)
        [average,log]=select_and_average(pooled,combined,groups{k},cfg);
        merged.(groups{k}).subjects{i}=average;
        merged.(groups{k}).trial_log{i}=log;
    end
    clear spectra pooled bhv
end
for k=1:numel(groups)
    merged.(groups{k})=summarize_subjects(merged.(groups{k}),subjectID,groups{k});
    for g=1:4
        name=cfg.condition_names{g};
        byCondition.(name).(groups{k})= ...
            summarize_subjects(byCondition.(name).(groups{k}),subjectID,groups{k});
    end
end
merged.pooling='Trial-weighted pooling before group-specific trimming.';
merged.valid_response='answer2 in {1,2} and finite rtime2 > 0';
merged.trim_rule='Separate within each trial group of each subject; merged and separate scopes independent.';
byCondition.valid_response=merged.valid_response;
byCondition.trim_rule=merged.trim_rule;
end

function [b,raw] = align_trials(s,raw,cfg,conditionIndex)
fields={'conflict','ischange','answer2','rtime2','oiti','osti2'};
assert(isfield(raw,'rejection_info'),'Missing rejection_info.');
n=numel(s.answer2);
for f=fields
    v=s.(f{1});
    assert((isnumeric(v)||islogical(v)) && isvector(v) && numel(v)==n,'Behavior field length mismatch.');
    s.(f{1})=v(:);
end
rejection=raw.rejection_info;
assert(iscell(rejection) && n==numel(rejection)*cfg.trials_per_run, ...
    'Run count and behavior length do not match trials_per_run.');
bad=[]; runID=[];
for r=1:numel(rejection)
    localBad=rejection{r}.bad_trials(:);
    assert(all(isfinite(localBad) & localBad>=1 & localBad<=cfg.trials_per_run & ...
        localBad==fix(localBad)),'Invalid run-local bad trial IDs.');
    localBad=unique(localBad);
    assert(rejection{r}.num_kept_trials==cfg.trials_per_run-numel(localBad), ...
        'num_kept_trials and bad_trials disagree.');
    bad=[bad;localBad+cfg.trials_per_run*(r-1)]; %#ok<AGROW>
    runID=[runID;repmat(r,rejection{r}.num_kept_trials,1)]; %#ok<AGROW>
end
kept=setdiff((1:n)',bad);
assert(numel(kept)==numel(raw.trial) && numel(raw.time)==numel(kept),'MEG/behavior trial count mismatch.');
b=struct;
for f=fields, b.(f{1})=s.(f{1})(kept); end
b.dur_iti=b.oiti-b.osti2+0.05;
b.original_trial=kept; b.run_idx=runID;
b.experimental_condition=repmat(conditionIndex,numel(kept),1);
% Explicit behavior columns replace opaque trialinfo in our computation copy.
raw.trialinfo=[b.conflict,b.ischange,b.answer2,b.original_trial,b.run_idx,b.experimental_condition];
raw=rmfield(raw,'rejection_info');
end

function freq = compute_spectrum(raw,cfg)
low=struct('output','pow','channel','MEGGRAD','method','mtmconvol', ...
    'taper','hanning','foi',cfg.low_foi,'toi',cfg.time_window(1):cfg.toi_step:cfg.time_window(2), ...
    't_ftimwin',repmat(cfg.low_t_ftimwin,1,numel(cfg.low_foi)), ...
    'keeptrials','yes','pad','nextpow2');
high=low; high.taper='dpss'; high.foi=cfg.high_foi;
high.t_ftimwin=repmat(cfg.high_t_ftimwin,1,numel(high.foi));
high.tapsmofrq=repmat(cfg.high_tapsmofrq,1,numel(high.foi));
a=ft_freqanalysis(low,raw); b=ft_freqanalysis(high,raw);
a=ft_combineplanar(struct('method','sum'),a);
b=ft_combineplanar(struct('method','sum'),b);
% Taper counts differ between Hanning and DPSS. Power is already estimated;
% do not concatenate taper bookkeeping as though it were a trial dimension.
for f={'cumtapcnt','cumsumcnt'}
    if isfield(a,f{1}), a=rmfield(a,f{1}); end
    if isfield(b,f{1}), b=rmfield(b,f{1}); end
end
freq=ft_appendfreq(struct('parameter','powspctrm','appenddim','freq'),a,b);
assert(strcmp(freq.dimord,'rpt_chan_freq_time'),'Unexpected spectral dimensions.');
% FieldTrip returns the closest FFT bins to requested frequencies.
assert(all(diff(freq.freq)>0) && round(freq.freq(1))==2 && ...
    round(freq.freq(end))==120,'Unexpected frequency grid.');
end

function [average,log] = select_and_average(freq,b,group,cfg)
valid=ismember(b.answer2,[1 2]) & isfinite(b.rtime2) & b.rtime2>0;
switch group
    case 'conflict', mask=valid & b.conflict==1;
    case 'congruent', mask=valid & b.conflict==0 & b.ischange==0;
    case 'left', mask=valid & b.answer2==1;
    case 'right', mask=valid & b.answer2==2;
end
candidate=find(mask);
removedLocal=[];
if cfg.trim_percent>0 && numel(candidate)>1
    removedLocal=extreme_indices(freq.powspctrm(candidate,:,:,:),cfg.trim_percent);
end
removed=candidate(removedLocal); retained=candidate;
retained(removedLocal)=[];
log=struct('n_selected',numel(candidate),'n_removed',numel(removed), ...
    'n_kept',numel(retained),'requested_percent',cfg.trim_percent, ...
    'actual_percent',0,'kept',[],'removed',[]);
if ~isempty(candidate), log.actual_percent=100*numel(removed)/numel(candidate); end
% Save original condition/run/trial IDs, even after pooling and trimming.
log.kept=table(b.experimental_condition(retained),b.run_idx(retained), ...
    b.original_trial(retained),b.answer2(retained),b.dur_iti(retained), ...
    'VariableNames',{'experimental_condition','run','original_trial','answer2','dur_iti'});
log.removed=table(b.experimental_condition(removed),b.run_idx(removed), ...
    b.original_trial(removed),'VariableNames',{'experimental_condition','run','original_trial'});
average=[];
if numel(retained)<cfg.min_trials
    warning('response_tfr:NoTrials','Group %s has %d retained trials; omitting this subject in this group.',group,numel(retained));
    return
end
trials=ft_selectdata(struct('trials',retained),freq);
average=ft_freqdescriptives(struct('keeptrials','no'),trials);
drop=intersect({'cfg','cumtapcnt','cumsumcnt','trialinfo'},fieldnames(average));
if ~isempty(drop), average=rmfield(average,drop); end
end

function removed = extreme_indices(power,percent)
% Same voxel-wise MAD / RMS scores as the original trim_by_percent_strict.
% 0% is handled by caller; positive percentages remove >=1, keeping >=1.
n=size(power,1); count=max(1,min(round(n*percent/100),n-1));
x=reshape(power,n,[]);
med=median(x,1,'omitnan'); mad=median(abs(x-med),1,'omitnan');
zero=mad==0; mad(zero)=std(x(:,zero),1,'omitnan')+eps;
z=abs(x-med)./(1.4826*mad);
scores=sqrt(mean(z.^2,2,'omitnan'));
% All-NaN trial spectra cannot be scored as ordinary low-power trials.
scores(isnan(scores))=Inf;
[~,order]=sort(scores,'descend');
removed=sort(order(1:count));
end

function item = summarize_subjects(item,subjectID,group)
present=~cellfun(@isempty,item.subjects);
item.valid_subject_mask=present;
item.all_subjects=subjectID;
item.subject=subjectID(present);
item.counts=table(string(subjectID), ...
    cellfun(@(x) x.n_selected,item.trial_log), ...
    cellfun(@(x) x.n_removed,item.trial_log), ...
    cellfun(@(x) x.n_kept,item.trial_log), ...
    'VariableNames',{'subject','selected','removed','kept'});
item.tfr=[];
if any(present)
    item.tfr=ft_freqgrandaverage(struct('keepindividual','yes'),item.subjects{present});
    item.tfr.subject=item.subject; item.tfr.condition=group;
end
item=rmfield(item,'subjects');
end
