%% decision cluster-permutation test
% cluster test，mutil test
cfg = [];
cfg.root_dir ='F:\xianliang\exp_data\OT_N2_MEG\MEG_source_408IndependDay\unitGain\decision_13_30Hz';
cfg.output_dir = fullfile(cfg.root_dir,'pairedF_cluster_results');
cfg.expected_subjects = 19;
cfg.expected_starts = -0.8:0.05:-0.1;
cfg.window_length = 0.2;
cfg.numrandomization = 5000;
cfg.clusteralpha = 0.05;
cfg.alpha = 0.05;
summary_table = run_sliding_source_F_cluster_to_brainnet(cfg);

% export NIfTI
raw_dir='F:\xianliang\exp_data\OT_N2_MEG\MEG_source_408IndependDay\unitGain\decision_13_30Hz\pairedF_cluster_results\statistics_mat';
stat_dir=dir([raw_dir,'\*.mat']);
output_nii='F:\xianliang\exp_data\OT_N2_MEG\MEG_source_408IndependDay\unitGain\decision_13_30Hz\pairedF_cluster_results\NIfTI_F_clusterSig1';
for i=1:length(stat_dir)
    % load(fullfile(raw_dir,stat_dir(i).name),'stat_res');
    export_cluster1_from_stat_file(fullfile(raw_dir,stat_dir(i).name),output_nii);
end




%% motor sensorimotor FDR test
% lr_path0='F:\xianliang\exp_data\OT_N2_MEG\MEG_source_408IndependDay\LeftRight_Source\F50-100Hz';
% lr_dir=dir(['F:\xianliang\exp_data\OT_N2_MEG\MEG_source_408IndependDay\LeftRight_Source','\*21.5hz(-0.2_0)*']);
% lr_path0='F:\xianliang\exp_data\OT_N2_MEG\MEG_source_408IndependDay\LeftRight_Source\F13-30Hz';
lr_path0='F:\xianliang\exp_data\OT_N2_MEG\MEG_source_408IndependDay\unitGain\left_right';
time_slide=-0.8:0.05:-0.1;
for i_slide=1:length(time_slide)
lr_dir=dir([lr_path0,'\*noDay*30-100*(',mat2str(time_slide(i_slide)),'_',mat2str(time_slide(i_slide)+0.2),')*']);
load(fullfile(lr_dir.folder,lr_dir.name,'all_source.mat'));

source_right_by_win{i_slide} = all_source_right;  % n_subject x 1 cell
source_left_by_win{i_slide}  = all_source_left;

end


win_start = -0.8:0.05:-0.1;
win_end   = win_start + 0.20;

cfg_fdr = [];
cfg_fdr.window_start = win_start;
cfg_fdr.window_end = win_end;

cfg_fdr.roi_method = 'aal';
cfg_fdr.ft_path = 'D:\software\matlab_toolbox\fieldtrip-20251218';
cfg_fdr.aal_include_postcentral = true;


cfg_fdr.ft_path = 'D:\software\matlab_toolbox\fieldtrip-20251218';
cfg_fdr.output_dir = 'F:\xianliang\exp_data\OT_N2_MEG\MEG_source_408IndependDay\unitGain\left_right\time_slide_niiData_all';
cfg_fdr.output_prefix = 'beta_motor_sliding_200ms_step50ms';

cfg_fdr.write_nifti = true;
cfg_fdr.write_interpolated_to_mri = true;

slide_res = sliding_motor_roi_fdr_from_source_sets( ...
    source_right_by_win, source_left_by_win, cfg_fdr);


%% motor gamma-band time slide source localization: sensorimotor F>2
raw_path='F:\xianliang\exp_data\OT_N2_MEG\MEG_source_408IndependDay\unitGain\left_right';
source_dir=dir([raw_path,'\*noDay*']);
for i=1:length(source_dir)
    load(fullfile(source_dir(i).folder,source_dir(i).name,'stat_cluster.mat'));
    cfg = [];

    cfg.ft_path = 'D:\software\matlab_toolbox\fieldtrip-20251218';

    cfg.roi_mode = 'sensorimotor';
    cfg.F_threshold = 2;
    cfg.connectivity = 26;
    mri_ref_file = fullfile(cfg.ft_path, 'template', 'anatomy', 'single_subj_T1.nii');

    output_dir =fullfile(source_dir(i).folder,source_dir(i).name);

    out = export_motor_roi_F_gt2_to_nifti(stat_res,mri_ref_file,output_dir,cfg);

end


