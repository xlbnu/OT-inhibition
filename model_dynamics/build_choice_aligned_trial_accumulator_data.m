function trialAccumulator = ...
        build_choice_aligned_trial_accumulator_data(trialBank,options)
%BUILD_CHOICE_ALIGNED_TRIAL_ACCUMULATOR_DATA GLM-ready trial trajectories.
% Observed ischange selects the matching model branch. Raw and pathwise
% baseline-corrected trajectories are retained for process101 and endpoint
% alignment, under fixed-B0 and dynamic-boundary normalization.

    arguments
        trialBank (1,1) struct
        options.minimumPaths (1,1) double = 30
        options.minimumBranchProbability (1,1) double = 0.01
        options.minimumReliableTrialsPerSubject (1,1) double = 1
    end
    assert(options.minimumPaths >= 1);
    assert(options.minimumBranchProbability >= 0 && ...
        options.minimumBranchProbability <= 1);

    nSubjects = size(trialBank.subject,1);
    nConditions = size(trialBank.subject,2);
    conditionRows = trialBank.meta.conditionRows;
    drugCode = [1 0 1 0];
    contextCode = [1 1 0 0];

    trialAccumulator.meta.created = datetime('now');
    trialAccumulator.meta.type = ...
        'trial-level choice-aligned accumulator data for dynamics GLM';
    trialAccumulator.meta.choiceDefinition = ...
        'observed ischange selects matching model stay/change branch';
    trialAccumulator.meta.chosenWinnerRelation = [ ...
        'chosen equals model winner only after conditioning on the ', ...
        'branch matching observed choice'];
    trialAccumulator.meta.subjectIDs = trialBank.meta.subjectIDs;
    trialAccumulator.meta.subjectRowsIn85 = ...
        trialBank.meta.subjectRowsIn85;
    trialAccumulator.meta.conditionRows = conditionRows;
    trialAccumulator.meta.conditionNames = trialBank.meta.conditionNames;
    trialAccumulator.meta.process101 = trialBank.meta.process101;
    trialAccumulator.meta.endStep = trialBank.meta.endStep;
    trialAccumulator.meta.dimensions = ...
        'each matrix is empirical trial x process point/model step';
    trialAccumulator.meta.branchProbabilityDefinition = ...
        'matching branch paths divided by all valid decision paths';
    trialAccumulator.meta.inputEvidenceScale = [ ...
        'Vself/Vother are confidence-derived model inputs before ', ...
        'subject-specific v0, k1, and k2 are applied'];
    trialAccumulator.meta.primaryGLMOutcome = ...
        'end.choiceAligned.baselineCorrected.B0.difference';
    trialAccumulator.config = options;
    trialAccumulator.subject = cell(nSubjects,nConditions);
    trialAccumulator.completed = false(nSubjects,nConditions);

    nEnd = numel(trialBank.meta.endStep);
    reliableTrialCounts = zeros(nSubjects,nConditions,nEnd);
    for c = 1:nConditions
        sourceCondition = conditionRows(c);
        for s = 1:nSubjects
            bank = trialBank.subject{s,c};
            out = initialize_trial_metadata_local(bank, ...
                drugCode(sourceCondition),contextCode(sourceCondition));
            out.process101.progress = bank.process101.progress;
            out.end.step = bank.end.step;
            out = copy_branch_data_local(out,bank);
            out = make_choice_aligned_local(out,bank,options);
            trialAccumulator.subject{s,c} = out;
            reliableTrialCounts(s,c,:) = reshape(sum( ...
                out.end.choiceAligned.isReliable,1),1,1,[]);
            trialAccumulator.completed(s,c) = true;
        end
    end

    trialAccumulator.support.end.nReliableTrials = ...
        uint16(reliableTrialCounts);
    trialAccumulator.support.end.nSubjectsWithAnyReliableTrial = ...
        uint16(squeeze(sum(reliableTrialCounts > 0,1)));
    trialAccumulator.support.end.nSubjectsMeetingMinimumTrials = ...
        uint16(squeeze(sum(reliableTrialCounts >= ...
        options.minimumReliableTrialsPerSubject,1)));
    trialAccumulator.meta.completed = datetime('now');
end

function out = initialize_trial_metadata_local(bank,drug,context)
    nTrials = numel(bank.observedChoice);
    choiceSign = ones(nTrials,1);
    choiceSign(bank.observedChoice == 2) = -1;
    vSelf = double(bank.stiSelf(:));
    vOther = double(bank.stiOther(:));
    vChosen = vSelf;
    vUnchosen = vOther;
    change = bank.observedChoice == 2;
    vChosen(change) = vOther(change);
    vUnchosen(change) = vSelf(change);
    confChosen = double(bank.confSelf1(:));
    confUnchosen = double(bank.confOther1(:));
    confChosen(change) = double(bank.confOther1(change));
    confUnchosen(change) = double(bank.confSelf1(change));
    pFinal = 0.5-1/32+double(bank.confSelf2(:))/16;
    pFinal = min(max(pFinal,eps),1-eps);

    out.trial.sourceTrialIndex = bank.sourceTrialIndex(:);
    out.trial.observedChoice = bank.observedChoice(:);
    out.trial.ischange = uint8(bank.observedChoice(:)-1);
    out.trial.choiceSignStayPositive = choiceSign;
    out.trial.drug = repmat(drug,nTrials,1);
    out.trial.context = repmat(context,nTrials,1);
    out.simulation.requestedPaths = uint32(bank.requestedPaths(:));
    out.simulation.validDecisionPaths = uint32(bank.validPaths(:));
    out.simulation.validDecisionProbability = single( ...
        double(bank.validPaths(:))./double(bank.requestedPaths(:)));
    out.simulation.nonDecisionProbability = single(1-double( ...
        out.simulation.validDecisionProbability));
    out.confidence.self1 = double(bank.confSelf1(:));
    out.confidence.other1 = double(bank.confOther1(:));
    out.confidence.self2 = double(bank.confSelf2(:));
    out.confidence.chosen1 = confChosen;
    out.confidence.unchosen1 = confUnchosen;
    out.confidence.finalProbability = pFinal;
    out.confidence.finalLogOdds = log(pFinal./(1-pFinal));
    out.input.Vself = vSelf;
    out.input.Vother = vOther;
    out.input.VsumSelfOther = vSelf+vOther;
    out.input.VdiffSelfOther = vSelf-vOther;
    out.input.Vchosen = vChosen;
    out.input.Vunchosen = vUnchosen;
    out.input.VsumChosen = vChosen+vUnchosen;
    out.input.VdiffChosen = vChosen-vUnchosen;
end

function out = copy_branch_data_local(out,bank)
    alignments = {'process101','end'};
    branches = {'stay','change'};
    normalizations = {'B0','boundary'};
    accumulators = {'self','other'};
    for a = 1:numel(alignments)
        alignment = alignments{a};
        for b = 1:numel(branches)
            branch = branches{b};
            out.(alignment).(branch).nPaths = ...
                bank.(alignment).(branch).nPaths;
            if strcmp(alignment,'process101')
                out.(alignment).(branch).branchProbability = ...
                    bank.process101.(branch).probability;
            end
            for n = 1:numel(normalizations)
                normalization = normalizations{n};
                for x = 1:numel(accumulators)
                    accumulator = accumulators{x};
                    out.(alignment).(branch).raw.(normalization). ...
                        (accumulator) = bank.(alignment).(branch). ...
                        (normalization).(accumulator);
                    out.(alignment).(branch).baselineCorrected. ...
                        (normalization).(accumulator) = ...
                        bank.(alignment).(branch).baselineCorrected. ...
                        (normalization).(accumulator);
                end
            end
        end
    end
end

function out = make_choice_aligned_local(out,bank,options)
    observedChoice = bank.observedChoice(:);
    nTrials = numel(observedChoice);
    stay = observedChoice == 1;
    change = observedChoice == 2;
    branchProbability = nan(nTrials,1,'single');
    branchProbability(stay) = bank.process101.stay.probability(stay);
    branchProbability(change) = ...
        bank.process101.change.probability(change);
    processPaths = zeros(nTrials,1);
    processPaths(stay) = bank.process101.stay.nPaths(stay);
    processPaths(change) = bank.process101.change.nPaths(change);
    endPaths = zeros(nTrials,size(bank.end.stay.nPaths,2));
    endPaths(stay,:) = bank.end.stay.nPaths(stay,:);
    endPaths(change,:) = bank.end.change.nPaths(change,:);
    out.process101.choiceAligned.nPaths = uint32(processPaths);
    out.process101.choiceAligned.branchProbability = branchProbability;
    out.process101.choiceAligned.isReliable = ...
        processPaths >= options.minimumPaths & ...
        branchProbability >= options.minimumBranchProbability;
    out.end.choiceAligned.nPaths = uint32(endPaths);
    out.end.choiceAligned.branchProbability = branchProbability;
    out.end.choiceAligned.isReliable = ...
        endPaths >= options.minimumPaths & ...
        branchProbability >= options.minimumBranchProbability;

    alignments = {'process101','end'};
    versions = {'raw','baselineCorrected'};
    normalizations = {'B0','boundary'};
    for a = 1:numel(alignments)
        alignment = alignments{a};
        for v = 1:numel(versions)
            version = versions{v};
            for n = 1:numel(normalizations)
                normalization = normalizations{n};
                staySelf = out.(alignment).stay.(version). ...
                    (normalization).self;
                stayOther = out.(alignment).stay.(version). ...
                    (normalization).other;
                changeSelf = out.(alignment).change.(version). ...
                    (normalization).self;
                changeOther = out.(alignment).change.(version). ...
                    (normalization).other;
                self = nan(size(staySelf),'single');
                other = nan(size(stayOther),'single');
                chosen = nan(size(staySelf),'single');
                unchosen = nan(size(stayOther),'single');
                self(stay,:) = staySelf(stay,:);
                other(stay,:) = stayOther(stay,:);
                chosen(stay,:) = staySelf(stay,:);
                unchosen(stay,:) = stayOther(stay,:);
                self(change,:) = changeSelf(change,:);
                other(change,:) = changeOther(change,:);
                chosen(change,:) = changeOther(change,:);
                unchosen(change,:) = changeSelf(change,:);
                target = struct();
                target.self = self;
                target.other = other;
                target.chosen = chosen;
                target.unchosen = unchosen;
                target.sum = chosen+unchosen;
                target.difference = chosen-unchosen;
                out.(alignment).choiceAligned.(version). ...
                    (normalization) = target;
            end
        end
    end
    out.process101.choiceAligned.nReliableTrials = uint16(sum( ...
        out.process101.choiceAligned.isReliable));
    out.end.choiceAligned.nReliableTrialsAtStep = uint16(sum( ...
        out.end.choiceAligned.isReliable,1));
end
