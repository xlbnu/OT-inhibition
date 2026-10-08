% === Script one: Anatomical modeling and saving  ===
base_sensor_path = 'F:\xianliang\exp_data\OT_N2_MEG\MEG-raw\OT_data';
base_mri_path='F:\xianliang\exp_data\OT_N2_MEG\MRI_NII';
% Load standard MNI 8mm grid
ft_path = 'D:\software\matlab_toolbox\fieldtrip-20251218';
template_grid = fullfile(ft_path, 'template', 'sourcemodel', 'standard_sourcemodel3d5mm.mat');
load(template_grid);
dir_mri = dir(base_mri_path);
dir_sensor=dir(base_sensor_path);dir_sensor(14)=[];
num_subjects = length(dir_mri)-2;
for sub = 1:num_subjects
    % 1. Read original T1 and Original MEG sensors
    sub_name = dir_mri(sub+2).name;sub_name2=regexprep(dir_sensor(sub+2).name,'-\d','');
    if ~isequal(sub_name,sub_name2)
        error('please check subject id');
    end
    dir_n1=dir(fullfile(base_mri_path,dir_mri(sub+2).name));
    is_dir=find(cat(1,dir_n1.isdir)==1);
    dir_n2=dir(fullfile(dir_n1(is_dir(3)).folder,dir_n1(is_dir(3)).name));
    mri_now=fullfile(dir_n2(3).folder,dir_n2(3).name);
    mri_orig    = ft_read_mri(mri_now);
    dir_s1=dir(fullfile(base_sensor_path,dir_sensor(sub+2).name));
    sensor_now=fullfile(dir_s1(3).folder,dir_s1(3).name,'n1_quat_trans_tsss.fif');
    sensor_grad = ft_convert_units(ft_read_sens(sensor_now), 'mm');
    
    % 2. Manually select landmarks (This opens a window)
    cfg_realign = []; cfg_realign.method = 'interactive'; cfg_realign.coordsys = 'neuromag';
    mri_aligned = ft_volumerealign(cfg_realign, mri_orig);
    

    % 3. Automatically compute the time-consuming head model, grid and lead field
    cfg_seg = []; cfg_seg.output = 'brain';
    mri_segmented = ft_volumesegment(cfg_seg, mri_aligned);
    
    cfg_headmodel = []; cfg_headmodel.method = 'singleshell';
    vol_individual = ft_prepare_headmodel(cfg_headmodel, mri_segmented);
    
    cfg_grid = []; cfg_grid.method = 'basedonmni'; cfg_grid.template = sourcemodel; 
    cfg_grid.nonlinear = 'yes'; cfg_grid.mri = mri_aligned; cfg_grid.unit = 'mm';
    grid_individual = ft_prepare_sourcemodel(cfg_grid);
  
    % figure check
    figure('Name', 'PERFECT INDIVIDUAL ALIGNMENT CHECK', 'Color', 'w', 'Position', [100 100 800 600]);
    hold on;
    % 1. Plot individual head model (Use the unit-converted vol_individual_mm)
    ft_plot_headmodel(vol_individual, 'facecolor', 'cortex', 'edgecolor', 'none', 'facealpha', 0.5);
    %
    % 2. Plot MEG sensors (Use the unit-converted sensor_grad_mm)
    ft_plot_sens(sensor_grad, 'style', 'ob', 'markersize', 4);
    % 3. Plot the warped grid's inside points (Use the unit-converted grid_individual_mm)
    if islogical(grid_individual.inside)
        inside_pos = grid_individual.pos(grid_individual.inside, :);
    else
        inside_pos = grid_individual.pos(grid_individual.inside > 0, :);
    end
    plot3(inside_pos(:,1), inside_pos(:,2), inside_pos(:,3), 'r.', 'MarkerSize', 6);
    % ft_plot_mesh(grid_individual.pos(grid_individual.inside,:));
    % Adjust the view
    view([-135 35]); axis equal; axis vis3d; axis off; camlight; lighting gouraud;
    title('PERFECT ALIGNMENT (AUTO MATCHED)', 'FontSize', 14, 'FontWeight', 'bold');
    hold off;
    input('请确认头模、源模、传感器位置对齐 (按回车接受结果):', 's');
    close(gcf); 

    % 4. Critical: Bundle all of this participant's"anatomical foundation"and save!
    % Note: Strongly recommended to add '-v7.3', because anatomical data can be large
    save_name = sprintf('F:\\xianliang\\exp_data\\OT_N2_MEG\\MRI_Preprocess\\%s_anatomy.mat', sub_name);
    save(save_name, 'mri_aligned', 'vol_individual', 'grid_individual', '-v7.3');
    
    fprintf('✔ 第 %d 号被试 %s 解剖模型构建并保存完毕！\n', sub,sub_name);
end

%%
base_anatomy_path = 'F:\xianliang\exp_data\OT_N2_MEG\MRI_Preprocess\head_model';
mri_dir=dir([base_anatomy_path,'\*.mat']);
base_sensor_path_ot = 'F:\xianliang\exp_data\OT_N2_MEG\MEG-raw\OT_data';
base_sensor_path_pl = 'F:\xianliang\exp_data\OT_N2_MEG\MEG-raw\PL_data';
dir_meg_ot=dir([base_sensor_path_ot,'\*-*']);
dir_meg_pl=dir([base_sensor_path_pl,'\*-*']);
sub_name_meg=cell(length(dir_meg_ot),2);
for i=1:length(dir_meg_ot)
    xx=strsplit(dir_meg_ot(i).name,'-');sub_name_meg{i,1}=xx{1};
end
for i=1:length(dir_meg_pl)
    xx=strsplit(dir_meg_pl(i).name,'-');sub_name_meg{i,2}=xx{1};
end
cond_dirs = {
'social-OT(250hz)','s',base_sensor_path_ot;
'social-PL(250hz)','s',base_sensor_path_pl;
'nonsocial-OT(250hz)','n',base_sensor_path_ot;
'nonsocial-PL(250hz)','n',base_sensor_path_pl;
};
base_lf_path='F:\xianliang\exp_data\OT_N2_MEG\MRI_Preprocess\lead_field_pressPutton_3run(250hz)';
if ~exist(base_lf_path, 'dir'), mkdir(base_lf_path); end
%
for i_sub=1:length(mri_dir)
    sub_name=regexprep(mri_dir(i_sub).name,'_anatomy.mat','');
    fprintf('\n=========================================\n');
    fprintf('▶ 正在处理被试: %s\n', sub_name);
    load(fullfile(mri_dir(i_sub).folder,mri_dir(i_sub).name));
    vol_individual = ft_convert_units(vol_individual, 'mm');
    grid_individual = ft_convert_units(grid_individual, 'mm');
    for n_cond=1:length(cond_dirs)
        cond_name=cond_dirs{n_cond,1};
        fprintf('▶ 正在处理condition: %s\n', cond_name);
        sub_idx=find(strcmp(sub_name_meg(:,1),sub_name));
        cond_name1=dir([cond_dirs{n_cond,3},['\',sub_name,'*']]);
        cond_name2=dir(fullfile(cond_name1.folder,cond_name1.name));
        leadfield_individual=cell(3,1);
        sensor_grad0=cell(3,1);
        for i_run=1:3
        cond_now=fullfile(cond_name2(3).folder,cond_name2(3).name,[cond_dirs{n_cond,2},num2str(i_run),'_quat_tsss.fif']);
        sensor_grad = ft_convert_units(ft_read_sens(cond_now), 'mm');
        xx_name=strsplit(cond_name1.name,'-');
        if ~isequal(xx_name{1},sub_name)
            error('please check subject id');
        end
        cfg_lf = []; cfg_lf.grad = sensor_grad; cfg_lf.headmodel = vol_individual;
        cfg_lf.sourcemodel = grid_individual; cfg_lf.channel = 'MEGGRAD';
        leadfield_individual_tmp = ft_prepare_leadfield(cfg_lf);

        leadfield_individual_tmp = rmfield(leadfield_individual_tmp, 'cfg');
        leadfield_individual{i_run,1}=leadfield_individual_tmp;
        sensor_grad0{i_run,1}=sensor_grad;
        end
        if ~exist(fullfile(base_lf_path,cond_name),'dir');mkdir(fullfile(base_lf_path,cond_name));end
        save_path = fullfile(base_lf_path,cond_name,sprintf('%s_leadfield_MEGGRAD.mat',sub_name));
        save(save_path, 'leadfield_individual','sensor_grad0', '-v7.3');
    end
end


%% ========================================================================
%  MEG Adaptive lead-field selector (Optimal Sensor Array Finder)
%  Function: Estimate multiple same-day runs' Run center displacement, Find the geometric center, and identify extreme outliers
% ========================================================================
base_sensor_path_ot = 'F:\xianliang\exp_data\OT_N2_MEG\MEG-raw\OT_data';
base_sensor_path_pl = 'F:\xianliang\exp_data\OT_N2_MEG\MEG-raw\PL_data';
dir_meg_ot=dir([base_sensor_path_ot,'\*-*']);
dir_meg_pl=dir([base_sensor_path_pl,'\*-*']);
dir_meg_ot(12)=[];dir_meg_pl(12)=[];
sub_name_meg=cell(length(dir_meg_ot),2);
for i=1:length(dir_meg_ot)
    xx=strsplit(dir_meg_ot(i).name,'-');sub_name_meg{i,1}=xx{1};
end
for i=1:length(dir_meg_pl)
    xx=strsplit(dir_meg_pl(i).name,'-');sub_name_meg{i,2}=xx{1};
end
if ~isequal(sub_name_meg(:,1),sub_name_meg(:,2)) | length(sub_name_meg)~=19
error('please check subject id');
end
name_head_pos={'s1_head_pos.pos','s2_head_pos.pos','s3_head_pos.pos';...
    'n1_head_pos.pos','n2_head_pos.pos','n3_head_pos.pos';};

for i_sub=1:length(dir_meg_ot)
    sub_name=sub_name_meg{i_sub,1};
    fprintf('\n=========================================\n');
    fprintf('▶ 正在处理被试: %s\n', sub_name);    
    ot_cond1=dir(fullfile(dir_meg_ot(i_sub).folder,dir_meg_ot(i_sub).name));
    pos_files_ot=cell(2,3);
    for i_pos=1:3
        for i_ctx=1:2
            pos_files_ot{i_ctx,i_pos}=fullfile(ot_cond1(3).folder,ot_cond1(3).name,name_head_pos{i_ctx,i_pos});
        end
    end
    pos_files_ot=pos_files_ot';
    ot_cond1=dir(fullfile(dir_meg_pl(i_sub).folder,dir_meg_pl(i_sub).name));
    pos_files_pl=cell(2,3);
    for i_pos=1:3
        for i_ctx=1:2
            pos_files_pl{i_ctx,i_pos}=fullfile(ot_cond1(3).folder,ot_cond1(3).name,name_head_pos{i_ctx,i_pos});
        end
    end
    pos_files_pl=pos_files_pl';
    %
    % [optimal_idx.os(i_sub,1), isolated_bad_run.os{i_sub,1},dist_matrix.os{i_sub,1}] = check_head_movement_mode_6run(pos_files_ot(:,1),'social-OT');
    % [optimal_idx.ps(i_sub,1), isolated_bad_run.ps{i_sub,1},dist_matrix.ps{i_sub,1}] = check_head_movement_mode_6run(pos_files_pl(:,1),'social-PL');
    % [optimal_idx.on(i_sub,1), isolated_bad_run.on{i_sub,1},dist_matrix.on{i_sub,1}] = check_head_movement_mode_6run(pos_files_ot(:,2),'nonsocial-OT');
    % [optimal_idx.pn(i_sub,1), isolated_bad_run.pn{i_sub,1},dist_matrix.pn{i_sub,1}] = check_head_movement_mode_6run(pos_files_pl(:,2),'nonsocial-PL');
    [optimal_idx.ot(i_sub,1), isolated_bad_run.ot{i_sub,1},dist_matrix.ot{i_sub,1}] = check_head_movement_mode_6run(pos_files_ot(:),'OT');
    [optimal_idx.pl(i_sub,1), isolated_bad_run.pl{i_sub,1},dist_matrix.pl{i_sub,1}] = check_head_movement_mode_6run(pos_files_pl(:),'PL');
end

save('F:\xianliang\exp_data\OT_N2_MEG\MRI_Preprocess\out_run_4condition\checkHead_idx_4condition.mat','optimal_idx','isolated_bad_run','dist_matrix');

function [optimal_idx, isolated_bad_runs, dist_matrix] = check_head_movement_mode_6run(pos_files, condition_str)
    % 🚨 Replace with all same-day runs for this participant Run 's .pos file paths
    num_runs = length(pos_files);
    mean_pos = zeros(num_runs, 3);
    valid_runs = true(num_runs, 1); % Track file validity
    
    fprintf('======================================================\n');
    fprintf('🔍 正在测算 [%s] %d 个 Run 的多维头位矩阵...\n', condition_str, num_runs);
    fprintf('======================================================\n');
    
    % 1. Extract all Run mean physical coordinates (X, Y, Z)
    for i = 1:num_runs
        if ~exist(pos_files{i}, 'file')
            warning('文件缺失: %s', pos_files{i});
            valid_runs(i) = false; continue;
        end
        data = importdata(pos_files{i});
        if isstruct(data), num_data = data.data; else, num_data = data; end
        % Convert to millimeters
        mean_pos(i, :) = mean(num_data(:, 5:7) * 1000, 1); 
    end
    
    % 2. Construct N x N Euclidean distance matrix (Distance Matrix)
    dist_matrix = zeros(num_runs, num_runs);
    for i = 1:num_runs
        for j = 1:num_runs
            if valid_runs(i) && valid_runs(j)
                dist_matrix(i, j) = norm(mean_pos(i,:) - mean_pos(j,:));
            end
        end
    end
    
    % --- Additional functionality: Print each Run pairwise displacement ---
    fprintf('\n📊 【各 Run 之间两两位移 (mm)】\n');
    for i = 1:num_runs
        for j = i+1:num_runs
            if valid_runs(i) && valid_runs(j)
                fprintf('  Run %d <-> Run %d: %6.2f mm\n', i, j, dist_matrix(i, j));
            end
        end
    end
    
    % 3. Find the geometric center (Medoid) —— The run with the shortest mean distance to all other Run runs
    mean_dist_to_others = zeros(num_runs, 1);
    for i = 1:num_runs
        if valid_runs(i)
            % Compute run i items Run distance to all other valid Run mean distance
            others = setdiff(find(valid_runs), i);
            if ~isempty(others)
                mean_dist_to_others(i) = mean(dist_matrix(i, others));
            else
                mean_dist_to_others(i) = inf;
            end
        else
            mean_dist_to_others(i) = inf; % Exclude invalid files
        end
    end
    [min_avg_dist, optimal_idx] = min(mean_dist_to_others);
    
    % --- New core logic: Find all runs whose distances to other Run are all > 10mm extreme outliers Run ---
    isolated_bad_runs = [];
    for i = 1:num_runs
        if valid_runs(i)
            % Get current Run distance to all other valid Run distances
            others = setdiff(find(valid_runs), i);
            if ~isempty(others)
                dists_to_others = dist_matrix(i, others);
                % If distances to all others exceed 10mm
                if all(dists_to_others > 10.0)
                    isolated_bad_runs = [isolated_bad_runs, i];
                end
            end
        end
    end
    
    % 4. Identify extreme outliers (Find maximum distance, For reporting only)
    max_dist_overall = max(dist_matrix(:));
    [worst_i, worst_j] = find(dist_matrix == max_dist_overall, 1);
    
    % ==================== Result output ====================
    fprintf('\n🏆 【最优导联场推荐】\n');
    fprintf('  最佳代表 Run: 第 %d 个 Run\n', optimal_idx);
    fprintf('  它与其他 Run 的平均偏差仅为: %.2f mm\n', min_avg_dist);
    fprintf('  👉 策略: 在 ft_appenddata 后，请使用 data_combined.grad = run_data{%d}.grad;\n', optimal_idx);
    
    fprintf('\n🚨 【极端位移排查】\n');
    fprintf('  全局最大位移: %.2f mm (发生在 Run %d 和 Run %d 之间)\n', max_dist_overall, worst_i, worst_j);
    
    if ~isempty(isolated_bad_runs)
        fprintf('\n❌ 警告: 发现与其他所有有效 Run 距离均超过 10mm 的孤立 Run！\n');
        fprintf('  诊断: 孤立坏 Run 的编号为 [ %s ]\n', num2str(isolated_bad_runs));
        fprintf('  👉 强烈建议: 在后续的 Append 和 Localizer 计算中，直接将这些 Run 剔除 (Drop)！\n');
        
        % Check whether enough remain after exclusion Run
        remaining_runs = sum(valid_runs) - length(isolated_bad_runs);
        if remaining_runs < 2
            fprintf('  ⚠️ 严重警告：剔除后剩余 Run 数量不足，无法进行稳定的合并计算，建议对该被试采用源空间独立平均法。\n');
        end
    elseif max_dist_overall > 10.0
        fprintf('\n⚠️ 警告: 虽然最大位移超过 10mm，但没有 Run 与**所有**其他 Run 完全脱节。\n');
        fprintf('  这通常意味着存在几个互相接近的小团体 (如 Run 1,2 接近，Run 4,5 接近，但它们之间隔得很远)。\n');
        fprintf('  👉 策略: 建议人工审查上面的位移矩阵，决定是拆分计算还是剔除某一部分。\n');
    else
        fprintf('\n✅ 状态完美，未发现导致模糊的极端位移，请放心合并全天数据。\n');
    end
    fprintf('======================================================\n');
end
