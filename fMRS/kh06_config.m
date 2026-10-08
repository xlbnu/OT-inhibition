function cfg = kh06_config()
% Configure paths for run_KH06.m.
% Edit these paths to point to the local, authorized study data.
package_dir = fileparts(mfilename('fullpath'));
cfg.data_dir = fullfile(package_dir, 'data');
cfg.behavior_file = fullfile(package_dir, 'data', 'behavior.mat');
end
