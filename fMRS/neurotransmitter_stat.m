

clear;

load('E:\xianliang\matlab_m\social_decision_m\data\OT-N2-MRS\result_data\OT_MRS_data.mat');
load('E:\xianliang\matlab_m\social_decision_m\data\OT-N2-MRS\sendout\Aggregated_MRS_Results.mat');
%
[gaba,glx,glu]=getNeurotransmitter_4Run(ots_data,pls_data,otn_data,pln_data,AllData);

codeRoot = fileparts(fileparts(mfilename('fullpath')));
outputDir = fullfile(codeRoot,'analysis_outputs','MRS_results');
if ~isfolder(outputDir), mkdir(outputDir); end
save(fullfile(outputDir,'Neurotransmitter_data.mat'),'gaba','glx'); % Aggregated concentrations

function [gaba, glx, glu] = getNeurotransmitter_4Run(ots_data, pls_data, otn_data, pln_data, AllData)
% Extract each participant's OT/PL, social/nonsocial two runs under each condition run concentrations.
% cf/cg: Second phase conflict/congruent; cf1/cg1: First phase.
% expOrder=1 Select the same-type run first two, expOrder=2 Select the last two.

data = {ots_data, pls_data, otn_data, pln_data};
prefixes = {'ots', 'pls', 'otn', 'pln'};
contextIndex = [1, 1, 2, 2];
metabolites = {'GABA', 'Glx', 'Glu'};
fitNames = {'Conflict_Phase2_Fit', 'Consistent_Phase2_Fit', ...
            'Conflict_Phase1_Fit', 'Consistent_Phase1_Fit'};
suffixes = {'cf', 'cg', 'cf1', 'cg1'};
nSubjects = numel(ots_data);

% Each output field is participant count × 2 runs.
outputs = {struct(), struct(), struct()};
for m = 1:3
    for d = 1:4
        for f = 1:4
            field = [prefixes{d}, '_', suffixes{f}];
            outputs{m}.(field) = nan(nSubjects, 2);
        end
    end
end

for i = 1:nSubjects
    % Retain the original participant-ID and experiment-order checks.
    for d = 1:2
        a = data{d}(i);
        b = data{d + 2}(i);
        if ~isequal(a.expOrder, b.expOrder) || ~isequal(a.mrs_num, b.mrs_num)
            error('please check mrs_number and drug_number');
        end
    end

    subjectName = sprintf('S%02d', ots_data(i).mrs_num);
    values = readConcentrations(AllData.(subjectName), fitNames, metabolites);

    for d = 1:4
        columns = (1:2) + 2 * (data{d}(i).expOrder - 1);
        for f = 1:4
            field = [prefixes{d}, '_', suffixes{f}];
            for m = 1:3
                outputs{m}.(field)(i, :) = reshape(values(columns, f, m, contextIndex(d)), 1, 2);
            end
        end
    end
end

[gaba, glx, glu] = outputs{:};
end

function values = readConcentrations(subjectData, fitNames, metabolites)
% Extract in original field order and Run1–Run4 order; NaN Retain in the corresponding run.
contexts = {'social', 'nonsocial'};
orderNames = fieldnames(subjectData);
values = nan(4, 4, 3, 2); % run × fit type × metabolite × context
runCount = [0, 0];

for order = 1:numel(orderNames)
    session = subjectData.(orderNames{order});
    for r = 1:4
        run = session.(sprintf('Run%d', r));
        names = fieldnames(run);
        if numel(names) ~= 1
            continue;
        end
        c = find(strcmp(contexts, names{1}), 1);
        if isempty(c)
            continue;
        end
        runCount(c) = runCount(c) + 1;
        for f = 1:4
            voxel = run.(contexts{c}).(fitNames{f}).vox1;
            for m = 1:3
                values(runCount(c), f, m, c) = voxel.(metabolites{m}).ConcIU_AlphaTissCorr;
            end
        end
    end
end
end
