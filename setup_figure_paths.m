function setup_figure_paths(fieldtripRoot,spmRoot)
% SETUP_FIGURE_PATHS Configure paths for manual figure reproduction.
%   From the package root: setup_figure_paths
%   For Fig.3/S6/S7: setup_figure_paths('/path/to/fieldtrip-20251218')
%   Brain maps in Fig.3/S5/S6/S7 also need SPM12:
%   setup_figure_paths('/path/to/fieldtrip-20251218','/path/to/spm12')
%   Alternatively, configure FieldTrip separately using its ft_defaults.
%   This function does not change the current directory or run analyses.
codeRoot = fileparts(mfilename('fullpath'));
addpath(fullfile(codeRoot,'utilities'));
addpath(fullfile(codeRoot,'figures'));
if nargin >= 2 && ~isempty(spmRoot)
    assert(isfile(fullfile(spmRoot,'spm.m')), ...
        'The supplied SPM root must contain spm.m.');
    addpath(spmRoot);
end
if nargin >= 1 && ~isempty(fieldtripRoot)
    assert(isfile(fullfile(fieldtripRoot,'ft_defaults.m')), ...
        'The supplied FieldTrip root must contain ft_defaults.m.');
    addpath(fieldtripRoot);
    ft_defaults;
elseif exist('ft_defaults','file') == 2
    ft_defaults;
end
fprintf('Figure scripts: %s\n',fullfile(codeRoot,'figures'));
fprintf('Bundled inputs: %s\n',fullfile(codeRoot,'figure_data'));
fprintf('Figure outputs: %s\n',fullfile(codeRoot,'figure_outputs'));
if exist('ft_freqgrandaverage','file') ~= 2
    fprintf('Fig.3/S6/S7 also require FieldTrip; supply its root as the argument.\n');
end
if exist('spm','file') ~= 2
    fprintf('Brain maps in Fig.3/S5/S6/S7 require SPM12; supply its root as the second argument.\n');
end
end
