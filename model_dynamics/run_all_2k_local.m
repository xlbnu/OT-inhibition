function run_all_2k_local()
%RUN_ALL_2K_LOCAL Validate inputs and build both two-k trajectory products.
% The original one-k workflow is not modified. This entry point runs the
% full-model EB workflow first, then the full-model GLM export, and finally
% the independently refitted no-lambda workflow.

    codeDir = fileparts(mfilename('fullpath'));
    addpath(codeDir);
    check_2k_trajectory_inputs();

    % The child scripts intentionally clear their own workspace, so resolve
    % the driver path afresh for every invocation instead of reusing a local
    % variable that a child script could remove.
    run(fullfile(fileparts(mfilename('fullpath')), ...
        'run_hierarchical_eb_85_to_meg19_2k.m'));
    run(fullfile(fileparts(mfilename('fullpath')), ...
        'run_export_trial_accumulator_glm_data_2k.m'));
    run(fullfile(fileparts(mfilename('fullpath')), ...
        'run_export_no_lambda_trial_accumulator_data_2k.m'));

    fprintf('All two-k local trajectory exports completed.\n');
end
