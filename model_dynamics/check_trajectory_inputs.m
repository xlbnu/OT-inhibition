function report = check_trajectory_inputs(dataDir)
%CHECK_TRAJECTORY_INPUTS Validate the current fit/data contract.
% This is a read-only check. It does not fit models or generate trajectories.

    if nargin < 1 || isempty(dataDir)
        codeRoot = fileparts(fileparts(mfilename('fullpath')));
        dataDir = fullfile(codeRoot,'figure_data');
    end
    dataDir = char(dataDir);
    behavior85 = fullfile(dataDir,'behavioral_results','bhv_data_85.mat');
    behavior19 = fullfile(dataDir,'behavioral_results','bhv_data_MEG_19.mat');
    full85 = fullfile(dataDir,'model_results','model_fit_parameter_85.mat');
    full19 = fullfile(dataDir,'model_results','model_fit_parameter_MEG19.mat');
    noLambda85 = fullfile(dataDir,'model_results','model_fit_parameter_85_NOLambda.mat');
    noLambda19 = fullfile(dataDir,'model_results','model_fit_parameter_MEG19_NOLambda.mat');

    files = {behavior85,behavior19,full85,full19,noLambda85,noLambda19};
    for i = 1:numel(files)
        assert(isfile(files{i}),'Missing input file: %s',files{i});
    end

    dataNames = {'ots_data','pls_data','otn_data','pln_data'};
    fullNames = {'ots_minb1','pls_minb1','otn_minb1','pln_minb1'};
    noLambdaNames = {'ots_minb2','pls_minb2','otn_minb2','pln_minb2'};

    d85 = load(behavior85,dataNames{:});
    d19 = load(behavior19,dataNames{:});
    pFull85 = load(full85,fullNames{:});
    pFull19 = load(full19,fullNames{:});
    pNoLambda85 = load(noLambda85,noLambdaNames{:});
    pNoLambda19 = load(noLambda19,noLambdaNames{:});

    subjectIDs85 = extract_ids_local(d85.ots_data);
    subjectIDs19 = extract_ids_local(d19.ots_data);
    [matched,rowsIn85] = ismember(subjectIDs19,subjectIDs85);
    assert(all(matched) && numel(unique(rowsIn85)) == 19, ...
        'MEG subject IDs do not map uniquely to the 85-subject cohort.');
    assert(isequal(rowsIn85(:)',38:56), ...
        'The MEG rows are not 38:56 in the current 85-subject order.');

    report = struct();
    report.dataDir = dataDir;
    report.subjectCount85 = numel(d85.ots_data);
    report.subjectCount19 = numel(d19.ots_data);
    report.megRowsIn85 = rowsIn85(:)';
    report.conditionNames = {'OT-social','PL-social', ...
        'OT-nonsocial','PL-nonsocial'};
    report.validConflictTrials85 = zeros(85,4);
    report.validConflictTrials19 = zeros(19,4);

    for c = 1:4
        d85c = d85.(dataNames{c});
        d19c = d19.(dataNames{c});
        f85 = pFull85.(fullNames{c});
        f19 = pFull19.(fullNames{c});
        n85 = pNoLambda85.(noLambdaNames{c});
        n19 = pNoLambda19.(noLambdaNames{c});

        assert(numel(d85c) == 85 && numel(f85) == 85 && ...
            numel(n85) == 85,'Condition %d is not 85 subjects.',c);
        assert(numel(d19c) == 19 && numel(f19) == 19 && ...
            numel(n19) == 19,'Condition %d is not 19 subjects.',c);

        for s = 1:85
            theta = get_parameter_value_local(f85(s));
            theta0 = get_parameter_value_local(n85(s));
            assert(numel(theta) == 9 && all(isfinite(theta)), ...
                'Full fit condition %d subject %d is not a finite 9-vector.',c,s);
            assert(numel(theta0) == 9 && all(isfinite(theta0)), ...
                'No-lambda fit condition %d subject %d is not a finite 9-vector.',c,s);
            assert(abs(theta0(7)) < 1e-12, ...
                'No-lambda fit condition %d subject %d has nonzero lambda.',c,s);
            report.validConflictTrials85(s,c) = count_valid_local( ...
                d85c(s),c <= 2);
        end
        for s = 1:19
            theta = get_parameter_value_local(f19(s));
            theta0 = get_parameter_value_local(n19(s));
            assert(numel(theta) == 9 && all(isfinite(theta)), ...
                'MEG full fit condition %d subject %d is not a finite 9-vector.',c,s);
            assert(numel(theta0) == 9 && all(isfinite(theta0)), ...
                'MEG no-lambda fit condition %d subject %d is not a finite 9-vector.',c,s);
            assert(abs(theta0(7)) < 1e-12, ...
                'MEG no-lambda fit condition %d subject %d has nonzero lambda.',c,s);
            report.validConflictTrials19(s,c) = count_valid_local( ...
                d19c(s),c <= 2);
        end
    end

    fprintf('input check passed.\n');
    fprintf('85 subjects: %d; MEG subjects: %d; MEG rows: %s\n', ...
        report.subjectCount85,report.subjectCount19,mat2str(rowsIn85(:)'));
    fprintf('Full model: 9 finite parameters in all conditions.\n');
    fprintf('No-lambda model: parameter 7 is zero in all conditions.\n');
    fprintf('MEG valid conflict trials per condition: %s\n', ...
        mat2str(sum(report.validConflictTrials19,1)));
end

function ids = extract_ids_local(data)
    ids = strings(numel(data),1);
    for s = 1:numel(data)
        raw = data(s).subsName;
        if iscell(raw) && isscalar(raw)
            raw = raw{1};
        end
        value = string(raw);
        ids(s) = value(1);
    end
end

function value = get_parameter_value_local(record)
    if isfield(record,'optParams_FullArray')
        value = double(record.optParams_FullArray);
    elseif isfield(record,'optParams')
        value = double(record.optParams);
    else
        error('No optParams or optParams_FullArray field found.');
    end
end

function count = count_valid_local(data,isSocial)
    if isSocial
        other = data.o1conf1(:);
    else
        other = data.inforate(:);
    end
    valid = data.conflict(:) == 1 & ...
        data.answer1(:) ~= -1 & data.answer2(:) ~= -1 & ...
        data.conf1(:) ~= -1 & data.conf2(:) ~= -1 & ...
        isfinite(data.conf1(:)) & isfinite(data.conf2(:)) & ...
        isfinite(other) & isfinite(data.ischange(:));
    count = sum(valid);
end
