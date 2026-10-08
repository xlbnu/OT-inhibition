% Behavioral analysis: accuracy/synergy, sigmoid curves, mixed-model betas.
% Run this script in MATLAB. Local functions are at the end of this file.
% All analyses, including pooled curves, use the same 85 retained subjects.
% Four-condition column order: OT-social, PL-social, OT-nonsocial, PL-nonsocial.
% AUC uses both initial and final confidence completeness, matching the original.
% rtime1 units and conf1 scale are kept unchanged.
% RT uses conflict==1 & answer2~=-1 for BOTH drugs. The old descriptive
% code applied this response criterion only to OT; PL omitted it.
% Consequently PL RT means can differ from the old saved results.
% Overall mixed-model betas include cg_so1 / cg_ns1 (congruent stay trials).
% Files are written to analysis_outputs/behavioral_results, preserving figure inputs.
% Overall outputs: osub, psub, synergy, sig_isc, sig_conf, sig_isc0,
% sig_conf0, beta_ot0, beta_pl0, choiceDeltaV, choiceSigmaV,
% confidenceDeltaV, confidenceSigmaV.
% Cohort outputs: bhv_ot0/bhv_pl0, MEG_ot0/MEG_pl0, MRS_ot0/MRS_pl0.
% sbeta_ot0/sbeta_pl0 concatenate the independently fitted cohort betas.
% beta_ot0/beta_pl0 are from the overall model, not the cohort fits.
% Curves: predict(sig_isc0{1},x) gives p(stay); confidence requires
% 1./(1+exp(-predict(sig_conf0{1},x))). Individual models use {subject,condition}.
% Unlike previously saved pooled models using excluded participants, here
% individual and pooled analyses use the same selected participant set.

%% Input all cohorts data
codeRoot = fileparts(fileparts(mfilename('fullpath')));
load(fullfile(codeRoot,'figure_data','behavioral_results','bhv_data_85.mat'));
cohortID = [ones(37,1);2*ones(19,1);3*ones(29,1)];
assert(all([numel(ots_data),numel(pls_data),numel(otn_data),numel(pln_data)]==85), 'Expected 85 aligned subjects.');
outputDir = fullfile(codeRoot,'analysis_outputs','behavioral_results');
if ~isfolder(outputDir), mkdir(outputDir); end

%% 1. Conflict-trial accuracy, synergy, AUC, RT, confidence and coherence
[osub,psub,synergy,coh_data] = getAccuracy(ots_data,pls_data,otn_data,pln_data);

save(fullfile(outputDir,'all_accauc.mat'),'osub','psub');
% Conflict-trial coherence is returned by getAccuracy.
save(fullfile(outputDir,'coherence_data.mat'),'coh_data');

% MEG 
[osub,psub,~,~] = getAccuracy(ots_data(38:56),pls_data(38:56),otn_data(38:56),pln_data(38:56));

save(fullfile(outputDir,'MEG_acc_auc.mat'),'osub','psub');


%% 2. Individual and pooled choice/confidence curves
[sig_isc,sig_conf,sig_isc0,sig_conf0] = getSigmoid(ots_data,pls_data,otn_data,pln_data);

save(fullfile(outputDir,'sigmoid_data.mat'),'sig_isc','sig_isc0','sig_conf','sig_conf0');

%% 3. Conflict and congruent mixed models and DeltaV / SigmaV coefficients
[beta_ot0,beta_pl0,effects] = getGLMbeta(ots_data,pls_data,otn_data,pln_data,cohortID);
choiceDeltaV = effects.choiceDeltaV;
choiceSigmaV = effects.choiceSigmaV;
confidenceDeltaV = effects.confidenceDeltaV;
confidenceSigmaV = effects.confidenceSigmaV;
congruentDeltaV = effects.congruentDeltaV;
congruentSigmaV = effects.congruentSigmaV;

save(fullfile(outputDir,'all_GLM_beta_withGroup_fixExpID.mat'),'beta_ot0','beta_pl0');

%% 4. cohort mixed models for each cohort
bhv_idx=1:37;
[bhv_ot0,bhv_pl0] = getGLMbeta(ots_data(bhv_idx),pls_data(bhv_idx),otn_data(bhv_idx),pln_data(bhv_idx),cohortID(bhv_idx));
MEG_idx=38:56;
[MEG_ot0,MEG_pl0] = getGLMbeta(ots_data(MEG_idx),pls_data(MEG_idx),otn_data(MEG_idx),pln_data(MEG_idx),cohortID(MEG_idx));
MRS_idx=57:85;
[MRS_ot0,MRS_pl0] = getGLMbeta(ots_data(MRS_idx),pls_data(MRS_idx),otn_data(MRS_idx),pln_data(MRS_idx),cohortID(MRS_idx));
ot_field=fieldnames(bhv_ot0);
for i=1:length(ot_field)
sbeta_ot0.(ot_field{i})=[bhv_ot0.(ot_field{i});MEG_ot0.(ot_field{i});MRS_ot0.(ot_field{i})];
end
pl_field=fieldnames(bhv_pl0);
for i=1:length(pl_field)
sbeta_pl0.(pl_field{i})=[bhv_pl0.(pl_field{i});MEG_pl0.(pl_field{i});MRS_pl0.(pl_field{i})];
end

save(fullfile(outputDir,'singleExp_GLM_beta.mat'),'sbeta_ot0','sbeta_pl0');
save(fullfile(outputDir,'MEG_GLM_beta.mat'),'MEG_ot0','MEG_pl0');
save(fullfile(outputDir,'MRS_GLM_beta.mat'),'MRS_ot0','MRS_pl0');




%% Local function 1: required conflict-trial descriptive statistics
function [osub,psub,synergy,coh_data] = getAccuracy(ots,pls,otn,pln)
groups={ots,pls,otn,pln}; n=numel(ots);
accSelf=cell(1,4); accOther=cell(1,4);
auc=cell(1,4); rt=cell(1,4); conf=cell(1,4); coh=cell(1,4);
for g=1:4
    assert(numel(groups{g})==n,'Condition subject counts differ.');
    accSelf{g}=nan(n,2); accOther{g}=nan(n,1);
    auc{g}=nan(n,2); rt{g}=nan(n,1); conf{g}=nan(n,1); coh{g}=nan(n,1);
    for i=1:n
        s=groups{g}(i); use=s.conflict(:)==1;
        correct1=s.isc1(:); correct2=s.isc2(:);
        conf1=s.conf1(:); conf2=s.conf2(:);
        % Original getDescriptiveStats AUC mask: initial correctness and
        % both confidence ratings must be nonmissing; use the same mask twice.
        aucMask=use & ~isnan(correct1) & ~isnan(conf1) & ~isnan(conf2);
        scorePair=[conf1,conf2]; correctPair=[correct1,correct2];
        for stage=1:2
            labels=correctPair(aucMask,stage)==1;
            scores=scorePair(aucMask,stage);
            npos=sum(labels); nneg=numel(labels)-npos;
            % Equivalent to rankauc; ties receive average ranks.
            if npos>0 && nneg>0
                ranks=tiedrank(scores);
                auc{g}(i,stage)=(sum(ranks(labels))-npos*(npos+1)/2)/(npos*nneg);
            end
        end
        responseTime=s.rtime1(:); finalAnswer=s.answer2(:); coherence=s.coh(:);
        % Apply the same response-validity rule in both drug conditions.
        rt{g}(i)=mean(responseTime(use & finalAnswer~=-1),'omitnan');
        conf{g}(i)=mean(conf1(use),'omitnan');
        coh{g}(i)=mean(coherence(use),'omitnan');
        accSelf{g}(i,1)=sum(correct1(use)==1)/sum(~isnan(correct1(use)));
        accSelf{g}(i,2)=sum(correct2(use)==1)/sum(~isnan(correct2(use)));
        if g<=2
            answer=s.o1answer1(:); target=s.dir(:);
            accOther{g}(i)=sum(answer(use)==target(use))/sum(~isnan(answer(use)));
        else
            correct=s.infoisc(:); answer=s.infoans(:);
            accOther{g}(i)=sum(correct(use)==1)/sum(~isnan(answer(use)));
        end
    end
end
osub=struct('acc_sc',accSelf{1},'acc_oc',accOther{1}, ...
    'acc_nc',accSelf{3},'acc_ic',accOther{3}, ...
    'auc_sc',auc{1},'auc_nc',auc{3},'rts_rdmcf',rt{1},'rtn_rdmcf',rt{3}, ...
    'conf_cfs',conf{1},'conf_cfn',conf{3});
psub=struct('acc_sc',accSelf{2},'acc_oc',accOther{2}, ...
    'acc_nc',accSelf{4},'acc_ic',accOther{4}, ...
    'auc_sc',auc{2},'auc_nc',auc{4},'rts_rdmcf',rt{2},'rtn_rdmcf',rt{4}, ...
    'conf_cfs',conf{2},'conf_cfn',conf{4});
coh_data=struct('os_cf',coh{1},'ps_cf',coh{2},'on_cf',coh{3},'pn_cf',coh{4});
synergy=nan(n,4);
for g=1:4
    synergy(:,g)=accSelf{g}(:,2)-max([accSelf{g}(:,1),accOther{g}],[],2);
end
end

%% Local function 2: individual and pooled curves
function [sig_isc,sig_conf,sig_isc0,sig_conf0] = getSigmoid(ots,pls,otn,pln)
% Choice predictor: self-minus-other log odds, response: stay.
% Confidence predictor: chosen-minus-unchosen, response: final log odds.
% Convert confidence predictions with 1./(1+exp(-predict(model,x))).
groups={ots,pls,otn,pln}; n=numel(ots);
sig_isc=cell(n,4); sig_conf=cell(n,4);
sig_isc0=cell(1,4); sig_conf0=cell(1,4);
for g=1:4
    assert(numel(groups{g})==n,'Condition subject counts differ.');
    pooled=cell(n,4);
    for i=1:n
        s=groups{g}(i); use=s.conflict(:)==1;
        x=s.lo_s1(:)-s.lo_o1(:); stay=s.ischange(:)==0;
        dc=s.lo_chosen(:)-s.lo_unchosen(:); y=s.lo_s2(:);
        x=x(use); stay=stay(use); dc=dc(use); y=y(use);
        sig_isc{i,g}=fitglm(x,stay,'Distribution','binomial','Link','logit');
        sig_conf{i,g}=fitglm(dc,y);
        pooled(i,:)={x,stay,dc,y};
    end
    sig_isc0{g}=fitglm(vertcat(pooled{:,1}),vertcat(pooled{:,2}), ...
        'Distribution','binomial','Link','logit');
    sig_conf0{g}=fitglm(vertcat(pooled{:,3}),vertcat(pooled{:,4}));
end
end

%% Local function 3: mixed-model subject coefficients
function [beta_ot0,beta_pl0,effects] = getGLMbeta(ots,pls,otn,pln,cohortID)
% Overall fit: cohort fixed covariate and cohort random intercept.
% Single-cohort fit: omit cohort terms (they are constant within a cohort).
% Cohort fits re-estimate models; they are not subsets of overall betas.
groups={ots,pls,otn,pln}; n=numel(ots);
if nargin<5 || isempty(cohortID), cohortID=ones(n,1); end
assert(isnumeric(cohortID) && numel(cohortID)==n && ...
    all(isfinite(cohortID(:))),'Invalid cohortID.');
cohortID=cohortID(:); multipleCohorts=numel(unique(cohortID))>1;
choice=cell(1,4); confidence=cell(1,4); congruent=cell(1,4);
for g=1:4
    assert(numel(groups{g})==n,'Condition subject counts differ.');
    fprintf('Mixed models: %d subjects, condition %d/4\n',n,g);
    for kind=1:3
        isChoice=kind==1;
        if isChoice
            predictors={'lo_plus','lo_minus'}; response='isstay';
        elseif kind==2
            predictors={'lo_cplus','lo_cminus'}; response='lo_s2';
        else
            predictors={'lo_plus','lo_minus'}; response='lo_s2';
        end
        pieces=cell(n,1);
        for i=1:n
            s=groups{g}(i);
            if kind==3
                use=s.conflict(:)==0 & s.ischange(:)==0;
            else
                use=s.conflict(:)==1;
            end
            if kind==1
                s.isstay=double(s.ischange==0);s.isstay(~isfinite(s.ischange(:))) = NaN;
            end
            y=s.(response)(:);
            x=[s.(predictors{1})(:),s.(predictors{2})(:)];
            y=y(use); x=x(use,:);
            if isChoice
                valid=~isnan(y) & ~any(isnan(x),2);
                y=y(valid); x=x(valid,:);
                x=zscore(x);
            else
                % Preserve original helper's response-first NaN removal.
                valid=~isnan(y);
                y=y(valid); x=x(valid,:);
                x=x-mean(x,1,'omitnan');
            end
            assert(~isempty(y),'Subject %d has no retained trials.',i);
            pieces{i}=table(y,repmat(i,numel(y),1),x(:,1),x(:,2), ...
                repmat(cohortID(i),numel(y),1), ...
                'VariableNames',{response,'id',predictors{1},predictors{2},'exp_id'});
        end
        tbl=vertcat(pieces{:});
        if multipleCohorts
            % Numeric exp_id matches the original overall analysis.
            formula=sprintf('%s~%s+%s+exp_id+(%s+%s|id)+(1|exp_id)', ...
                response,predictors{1},predictors{2},predictors{1},predictors{2});
        else
            formula=sprintf('%s~%s+%s+(%s+%s|id)', ...
                response,predictors{1},predictors{2},predictors{1},predictors{2});
        end
        if isChoice
            model=fitglme(tbl,formula,'Distribution','Binomial','Link','logit');
        else
            model=fitlme(tbl,formula);
        end
        [re,~,stats]=randomEffects(model); fixed=fixedEffects(model);
        b=nan(n,3); terms={'(Intercept)',predictors{1},predictors{2}};
        for i=1:n
            for j=1:3
                ri=strcmp(stats.Group,'id') & str2double(string(stats.Level))==i & ...
                    strcmp(stats.Name,terms{j});
                fi=strcmp(model.CoefficientNames,terms{j});
                assert(sum(ri)==1 && sum(fi)==1,'Cannot map subject coefficient.');
                b(i,j)=fixed(fi)+re(ri);
            end
        end
        if isChoice
            choice{g}=b;
        elseif kind==2
            confidence{g}=b;
        else
            congruent{g}=b;
        end
    end
end
beta_ot0=struct('isc_so1',choice{1},'isc_ns1',choice{3}, ...
    'cf_so20',confidence{1},'cf_ns20',confidence{3}, ...
    'cg_so1',congruent{1},'cg_ns1',congruent{3});
beta_pl0=struct('isc_so1',choice{2},'isc_ns1',choice{4}, ...
    'cf_so20',confidence{2},'cf_ns20',confidence{4}, ...
    'cg_so1',congruent{2},'cg_ns1',congruent{4});
effects=struct;
effects.columnOrder={'OT_social','PL_social','OT_nonsocial','PL_nonsocial'};
effects.choiceDeltaV=nan(n,4); effects.choiceSigmaV=nan(n,4);
effects.confidenceDeltaV=nan(n,4); effects.confidenceSigmaV=nan(n,4);
effects.congruentDeltaV=nan(n,4); effects.congruentSigmaV=nan(n,4);
for g=1:4
    effects.choiceDeltaV(:,g)=choice{g}(:,3);
    effects.choiceSigmaV(:,g)=choice{g}(:,2);
    effects.confidenceDeltaV(:,g)=confidence{g}(:,3);
    effects.confidenceSigmaV(:,g)=confidence{g}(:,2);
    effects.congruentDeltaV(:,g)=congruent{g}(:,3);
    effects.congruentSigmaV(:,g)=congruent{g}(:,2);
end
end
