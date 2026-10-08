function run_all_local()
%RUN_ALL_LOCAL Validate inputs and build both trajectory products.
% This entry point runs the
% full-model EB workflow first, then the full-model GLM export, and finally
% the independently refitted no-lambda workflow.

    codeDir = fileparts(mfilename('fullpath'));
    addpath(codeDir);
    check_trajectory_inputs();

    % The child scripts intentionally clear their own workspace, so resolve
    % the driver path afresh for every invocation instead of reusing a local
    % variable that a child script could remove.
    run(fullfile(fileparts(mfilename('fullpath')), ...
        'run_hierarchical_eb_85_to_meg19.m'));
    run(fullfile(fileparts(mfilename('fullpath')), ...
        'run_export_trial_accumulator_glm_data.m'));
    run(fullfile(fileparts(mfilename('fullpath')), ...
        'run_export_no_lambda_trial_accumulator_data.m'));

    fprintf('All local trajectory exports completed.\n');
end
