%% Direct 4 items condition comparison
% Standardized equivalence bound for comparisons among conditions within each parameter.
sd_multiplier = 0.6;

%% full model
[pair_table, summary_table] = equivalence_four_conditions_sd( ...
    ots_xx1, pls_xx1, otn_xx1, pln_xx1, ...
    'SDMultiplier', sd_multiplier, ...
    'ParameterOrder', [7 1:6 8 9], ...
    'Alpha', 0.05);

%% two lambda model
[pair_table, summary_table] = equivalence_four_conditions_sd( ...
    ots_xx1, pls_xx1, otn_xx1, pln_xx1, ...
    'SDMultiplier', sd_multiplier, ...
    'ParameterOrder', [7 8 1:6 9 10], ...
    'Alpha', 0.05);





%%  equivalence for parameter levels.
% Load aligned subjects x parameters matrices for ONE model before running:
% ots_xx1, pls_xx1, otn_xx1, pln_xx1.

model_n_parameters = 10;  % Use 9 for the symmetric model or 10 for two-direction inhibition.
required_conditions = 4;
sd_multiplier = 0.75;   % Standardized bound for comparisons between different parameters.
alpha = 0.05;

switch model_n_parameters
    case 9
        pn = [7 1:6 8 9];
        parameter_min = [1e-5, 0.001, 0, 0, 0.001,0.001, 0.1, 10, 10];
        parameter_max = [0.15 1 1 1 10 10 20 1000 1000];
        parameter_names = {'lambda','v0','bs/b0','bo/b0', ...
            'ks','kp','sigma','tau','t0'};
    case 10
        pn = [7 8 1:6 9 10];
        parameter_min = [1e-5, 1e-5,0.001, 0, 0, 0.001,0.001, 0.1, 10, 10];
        parameter_max = [0.15 0.15 1 1 1 10 10 20 1000 1000];
        parameter_names = {'lambda_o_to_s','lambda_s_to_o', ...
            'v0','bs/b0','bo/b0','ks','kp','sigma','tau','t0'};
    otherwise
        error('Choose model_n_parameters = 9 or 10.');
end
assert(size(ots_xx1,2)==model_n_parameters, ...
    'Load the correct model data before selecting model_n_parameters.');

[condition_table, partial_summary, partial_details] = ...
    equivalence_parameter_levels_partial_conjunction_sd( ...
    ots_xx1,pls_xx1,otn_xx1,pln_xx1, ...
    'RequiredConditions',required_conditions, ...
    'ParameterOrder',pn,'ParameterMin',parameter_min, ...
    'ParameterMax',parameter_max,'ParameterNames',parameter_names, ...
    'SDMultiplier',sd_multiplier,'Alpha',alpha);

% Primary 3/4 result after Holm correction across parameter pairs.
disp(partial_summary(:,{ ...
    'ParameterName1','ParameterName2','N', ...
    'P_TOST_OTS','P_TOST_PLS','P_TOST_OTN','P_TOST_PLN', ...
    'P_UthOrdered','P_PartialConjunction', ...
    'P_PartialConjunction_Holm','EquivalentAtLeastU_Holm', ...
    'P_All4_Strict_Holm','EquivalentAll4_Strict_Holm'}));

equivalent_at_least_3 = partial_summary( ...
    partial_summary.EquivalentAtLeastU_Holm==1,:);
disp(equivalent_at_least_3);

% Compare outcomes: supported for at least 3/4, but not strict 4/4.
three_not_four = partial_summary( ...
    partial_summary.EquivalentAtLeastU_Holm==1 & ...
    partial_summary.EquivalentAll4_Strict_Holm==0,:);
disp(three_not_four(:,{ ...
    'ParameterName1','ParameterName2','LargestPCondition', ...
    'LargestConditionP','P_PartialConjunction_Holm', ...
    'P_All4_Strict_Holm'}));

% Matrices use pn / parameter_names order; diagonal is NaN.
p_partial_holm_matrix = partial_details.P_PartialConjunction_Holm_Matrix;
equivalent_3of4_matrix = partial_details.EquivalentAtLeastU_Holm_Matrix;

figure;imagesc(p_partial_holm_matrix);















function [condition_table, summary_table, details] = ...
    equivalence_parameter_levels_partial_conjunction_sd(ots, pls, otn, pln, varargin)
%EQUIVALENCE_PARAMETER_LEVELS_PARTIAL_CONJUNCTION_SD At least u/4 conditions.
%
% For every parameter pair, this function first performs four standardized
% paired TOSTs, one in each condition (OTS, PLS, OTN, PLN), by calling
% equivalence_parameter_levels_sd. It then tests whether equivalence is
% supported in at least u of the four conditions using a conservative
% Bonferroni partial-conjunction p value:
%
%   sort p_(1) <= ... <= p_(4)
%   P_PC(u/4) = min(1, (4-u+1) * p_(u)).
%
% For u=3: P_PC(3/4) = min(1, 2*p_(3)). The smallest u-1 p values are not
% sufficient by themselves. Rejecting the partial-conjunction null supports
% equivalence in AT LEAST u conditions, but does not identify which ones.
%
% Holm correction is applied to the P_PC values across all parameter pairs
% selected in this function call (28 pairs for 8 parameters; 36 for 9).
% The strict all-four result is also retained as a sensitivity comparison:
% P_all4 = max(the four condition TOST p values), followed by Holm across
% parameter pairs.
%
% Required inputs follow equivalence_parameter_levels_sd:
%   ParameterOrder : analysis-order column indices
%   ParameterMin   : fixed lower bounds in ParameterOrder
%   ParameterMax   : fixed upper bounds in ParameterOrder
%   ParameterNames : names in ParameterOrder
%   SDMultiplier   : symmetric standardized equivalence bound c
%   Alpha          : default .05
%
% Assumptions and limitations:
% - independent participants and normally distributed paired differences;
% - fixed bounds must be chosen independently of these results;
% - all four condition tests must be computable for a parameter pair;
% - a 3/4 claim permits one condition to be non-equivalent or inconclusive;
% - condition-specific p values are unadjusted component results.

ip = inputParser;
addParameter(ip,'RequiredConditions',3,@(x) isnumeric(x)&&isscalar(x)&& ...
    isfinite(x)&&x==fix(x)&&x>=1&&x<=4);
addParameter(ip,'SDMultiplier',[],@positiveScalar);
addParameter(ip,'Alpha',0.05,@(x) positiveScalar(x)&&x<0.5);
addParameter(ip,'ParameterOrder',[],@(x) isnumeric(x)&&isreal(x)&&isvector(x)&& ...
    numel(x)>=2&&all(isfinite(x))&&all(x>=1)&&all(x==fix(x)));
addParameter(ip,'ParameterMin',[],@(x) isnumeric(x)&&isreal(x)&&isvector(x)&& ...
    all(isfinite(x)));
addParameter(ip,'ParameterMax',[],@(x) isnumeric(x)&&isreal(x)&&isvector(x)&& ...
    all(isfinite(x)));
addParameter(ip,'ParameterNames',{},@(x) iscellstr(x)||isstring(x));
parse(ip,varargin{:});

u = ip.Results.RequiredConditions;
c = ip.Results.SDMultiplier;
alpha = ip.Results.Alpha;
pn = ip.Results.ParameterOrder(:)';
minima = double(ip.Results.ParameterMin(:)');
maxima = double(ip.Results.ParameterMax(:)');
names = cellstr(string(ip.Results.ParameterNames));
assert(~isempty(c)&&~isempty(pn)&&~isempty(minima)&&~isempty(maxima), ...
    'Specify SDMultiplier, ParameterOrder, ParameterMin, and ParameterMax.');

% Obtain valid condition-specific TOST p values and the strict 4/4 result.
[condition_table, strict_table, core] = equivalence_parameter_levels_sd( ...
    ots,pls,otn,pln,'SDMultiplier',c,'Alpha',alpha, ...
    'ParameterOrder',pn,'ParameterMin',minima,'ParameterMax',maxima, ...
    'ParameterNames',names);

pCondition = core.ConditionP;
[m,nCondition] = size(pCondition);
assert(nCondition==4,'This procedure requires exactly four conditions.');
complete = all(isfinite(pCondition),2);
pSorted = sort(pCondition,2,'ascend');
pUth = nan(m,1);
pPC = nan(m,1);
pUth(complete) = pSorted(complete,u);
pPC(complete) = min(1,(nCondition-u+1).*pUth(complete));

% Descriptive only: number of component TOST p values below nominal alpha.
nPassing = nan(m,1);
nPassing(complete) = sum(pCondition(complete,:)<alpha,2);

% Holm adjustment across all planned parameter-pair partial-conjunction tests.
pPCHolm = holmPlanned(pPC,complete);
eqPCRaw = nan(m,1); eqPCHolm = nan(m,1);
eqPCRaw(complete) = double(pPC(complete)<alpha);
eqPCHolm(complete) = double(pPCHolm(complete)<alpha);

% Diagnose the largest condition p value; it is not the 3/4 test statistic.
[largestP,largestIndex] = max(pCondition,[],2);
conditionNames = string(core.ConditionNames(:));
largestCondition = strings(m,1);
largestCondition(complete) = conditionNames(largestIndex(complete));
largestCondition(~complete) = missing;

% Preserve the original pair ordering from the strict four-condition table.
summary_table = table( ...
    strict_table.Parameter1,strict_table.Parameter2, ...
    strict_table.ParameterName1,strict_table.ParameterName2,strict_table.N, ...
    repmat(u,m,1),repmat(nCondition,m,1),complete,nPassing, ...
    pCondition(:,1),pCondition(:,2),pCondition(:,3),pCondition(:,4), ...
    pUth,pPC,eqPCRaw,pPCHolm,eqPCHolm,largestCondition,largestP, ...
    strict_table.GlobalP,strict_table.GlobalPHolm,strict_table.EquivalentHolm, ...
    'VariableNames',{ ...
    'Parameter1','Parameter2','ParameterName1','ParameterName2','N', ...
    'RequiredConditions','TotalConditions','AllConditionsComputed', ...
    'NConditionsPassingNominal','P_TOST_OTS','P_TOST_PLS','P_TOST_OTN', ...
    'P_TOST_PLN','P_UthOrdered','P_PartialConjunction', ...
    'EquivalentAtLeastU_Unadjusted','P_PartialConjunction_Holm', ...
    'EquivalentAtLeastU_Holm','LargestPCondition','LargestConditionP', ...
    'P_All4_Strict','P_All4_Strict_Holm','EquivalentAll4_Strict_Holm'});

% Symmetric matrices follow ParameterOrder; diagonals were not tested.
np = numel(pn);
pPCMatrix = nan(np); pPCHolmMatrix = nan(np); eqPCMatrix = nan(np);
pairs = core.PairIndices;
for h=1:m
    a=pairs(h,1); b=pairs(h,2);
    pPCMatrix(a,b)=pPC(h); pPCMatrix(b,a)=pPC(h);
    pPCHolmMatrix(a,b)=pPCHolm(h); pPCHolmMatrix(b,a)=pPCHolm(h);
    eqPCMatrix(a,b)=eqPCHolm(h); eqPCMatrix(b,a)=eqPCHolm(h);
end

metadata = struct('Alpha',alpha,'DzBounds',[-c c], ...
    'RequiredConditions',u,'TotalConditions',nCondition, ...
    'PartialConjunctionMethod','Bonferroni partial conjunction', ...
    'Formula',sprintf('min(1,(4-%d+1)*p_(%d))',u,u), ...
    'HolmFamilySize',m,'ParameterOrder',pn,'ParameterMin',minima, ...
    'ParameterMax',maxima,'ParameterNames',{names}, ...
    'Interpretation',sprintf('Equivalence in at least %d of 4 conditions',u), ...
    'DoesNotIdentifyConditions',true);
condition_table.Properties.UserData = metadata;
summary_table.Properties.UserData = metadata;
details = metadata;
details.ConditionNames = core.ConditionNames;
details.ConditionP = pCondition;
details.SortedConditionP = pSorted;
details.PairIndices = pairs;
details.ValidSubjectsByPair = core.ValidSubjectsByPair;
details.P_PartialConjunction_Matrix = pPCMatrix;
details.P_PartialConjunction_Holm_Matrix = pPCHolmMatrix;
details.EquivalentAtLeastU_Holm_Matrix = eqPCMatrix;
details.StrictAll4 = strict_table;

fprintf(['Partial conjunction: equivalence in at least %d/4 conditions; ', ...
    'P_PC=min(1,%d*p_(%d)).\n'],u,nCondition-u+1,u);
fprintf('Holm adjustment across %d parameter-pair P_PC values.\n',m);
end

function adjusted = holmPlanned(p,complete)
% Retain uncomputable members in the planned family with internal p=1.
m = numel(p);
work = p(:);
work(~isfinite(work)) = 1;
[sorted,index] = sort(work);
adjustedSorted = min(1,cummax((m:-1:1)'.*sorted));
adjusted = nan(m,1);
adjusted(index) = adjustedSorted;
adjusted(~complete) = NaN;
end

function ok=positiveScalar(x)
ok=isnumeric(x)&&isreal(x)&&isscalar(x)&&isfinite(x)&&x>0;
end




function [condition_table, summary_table, details] = equivalence_parameter_levels_sd(ots, pls, otn, pln, varargin)
%EQUIVALENCE_PARAMETER_LEVELS_SD Parameter-level equivalence in ALL conditions.
% Inputs: aligned subjects x parameters matrices, not stacked observations.
% Scale each parameter as (X - ParameterMin)/(ParameterMax - ParameterMin).
% Set minima to 10 for tau/t0 and 0 for all other parameters.
% ParameterMin must be supplied explicitly. For each parameter pair, test
% theta = population mean(D)/population SD(D), D = scaled P - scaled Q,
% separately in OTS, PLS, OTN, PLN. GlobalP = max(four TOST p-values).
% Holm correction is across all parameter-pair GlobalP values in this call.
% Condition p-values are unadjusted components, not separate adjusted claims.
% Requires independent subjects, normal paired differences, and Statistics
% and Machine Learning Toolbox. SDMultiplier MUST be explicitly specified.
% ParameterMin, ParameterMax and ParameterNames follow ParameterOrder.
% Complete cases are selected per parameter pair across ALL four conditions.
% Ordinary 90%% mean-difference CIs are descriptive, not CIs for theta.
% Flags: 1 = established, 0 = not established, NaN = unavailable.

ip = inputParser;
addParameter(ip,'SDMultiplier',[],@positiveScalar);
addParameter(ip,'Alpha',0.05,@(x) positiveScalar(x)&&x<0.5);
addParameter(ip,'ParameterOrder',[],@(x) isnumeric(x)&&isreal(x)&&isvector(x)&& ...
    numel(x)>=2&&all(isfinite(x))&&all(x>=1)&&all(x==fix(x)));
addParameter(ip,'ParameterMax',[],@(x) isnumeric(x)&&isreal(x)&&isvector(x)&& ...
    all(isfinite(x))&&all(x>0));
addParameter(ip,'ParameterMin',[],@(x) isnumeric(x)&&isreal(x)&&isvector(x)&& ...
    all(isfinite(x)));
addParameter(ip,'ParameterNames',{},@(x) iscellstr(x)||isstring(x));
parse(ip,varargin{:});
c = ip.Results.SDMultiplier;
alpha = ip.Results.Alpha;
pn = ip.Results.ParameterOrder(:)';
maxima = double(ip.Results.ParameterMax(:)');
minima = double(ip.Results.ParameterMin(:)');
assert(~isempty(c)&&~isempty(pn)&&~isempty(maxima)&&~isempty(minima), ...
    'Specify SDMultiplier, ParameterOrder, ParameterMin and ParameterMax explicitly.');
inputs = {ots,pls,otn,pln};
for j=1:4
    validateattributes(inputs{j},{'numeric'},{'real','2d'});
end
assert(isequal(size(ots),size(pls),size(otn),size(pln)), ...
    'All input matrices must have identical sizes and aligned subjects.');
assert(max(pn)<=size(ots,2)&&numel(unique(pn))==numel(pn), ...
    'ParameterOrder must contain unique existing column indices.');
np = numel(pn);
assert(numel(maxima)==np,'ParameterMax must follow ParameterOrder.');
assert(numel(minima)==np,'ParameterMin must follow ParameterOrder.');
assert(all(maxima>minima),'Each ParameterMax must exceed ParameterMin.');
labels = cellstr(string(ip.Results.ParameterNames));
if isempty(labels)
    labels = arrayfun(@(x) sprintf('Column%d',x),pn,'UniformOutput',false);
end
assert(numel(labels)==np,'ParameterNames must follow ParameterOrder.');
conditions = {'OTS','PLS','OTN','PLN'};
scaled = cell(1,4);
for j=1:4
    scaled{j} = bsxfun(@rdivide, ...
        bsxfun(@minus,double(inputs{j}(:,pn)),minima),maxima-minima);
end
pairs = nchoosek(1:np,2);
m = size(pairs,1);
rows = cell(4*m,20);
globalP = nan(m,1); counts = zeros(m,1); allComputed = false(m,1);
conditionP = nan(m,4); validSubjects = cell(m,1);

for h=1:m
    a=pairs(h,1); b=pairs(h,2);
    valid=true(size(ots,1),1);
    for j=1:4
        valid=valid & all(isfinite(scaled{j}(:,[a b])),2);
    end
    validSubjects{h}=valid;
    n=sum(valid); counts(h)=n;
    for j=1:4
        mu=NaN; sd=NaN; se=NaN; dz=NaN; t=NaN;
        lo=NaN; hi=NaN; pL=NaN; pU=NaN; pt=NaN; eq=NaN;
        status='Too few complete participants';
        if n>=2
            d=scaled{j}(valid,a)-scaled{j}(valid,b);
            mu=mean(d); sd=std(d,0);
            status='Zero or nonfinite SD: standardized test unavailable';
            if isfinite(mu)&&isfinite(sd)&&sd>0
                se=sd/sqrt(n); dz=mu/sd; t=sqrt(n)*dz;
                pL=nctcdf(t,n-1,-c*sqrt(n),'upper');
                pU=nctcdf(t,n-1,c*sqrt(n));
                hw=tinv(1-alpha,n-1)*se;
                lo=mu-hw; hi=mu+hw;
                if all(isfinite([pL pU]))
                    pt=max(pL,pU); eq=double(pt<alpha); status='Computed';
                else
                    status='Noncentral t calculation unavailable';
                end
            end
        end
        conditionP(h,j)=pt;
        rows((h-1)*4+j,:)={pn(a),pn(b),labels{a},labels{b},conditions{j}, ...
            n,n-1,mu,sd,se,dz,t,lo,hi,-c,c,pL,pU,pt,eq};
        statuses{(h-1)*4+j,1}=status; %#ok<AGROW>
    end
    allComputed(h)=all(isfinite(conditionP(h,:)));
    if allComputed(h)
        globalP(h)=max(conditionP(h,:));
    end
end
condition_table=cell2table(rows,'VariableNames', ...
    {'Parameter1','Parameter2','ParameterName1','ParameterName2','Condition', ...
    'N','DF','MeanDifference','SD_Difference','SE_Difference','Dz','T_Observed', ...
    'MeanCI_Lower_Unadjusted','MeanCI_Upper_Unadjusted','DzLowerBound', ...
    'DzUpperBound','P_Lower','P_Upper','P_TOST_Dz','EquivalentUnadjusted'});
condition_table.Status=statuses;

% Keep unavailable comparisons in the planned Holm family using internal p=1.
workP=globalP; workP(~isfinite(workP))=1;
[sortedP,ix]=sort(workP);
adjusted=min(1,cummax((m:-1:1)'.*sortedP));
holmP=nan(m,1); holmP(ix)=adjusted; holmP(~allComputed)=NaN;
eqRaw=nan(m,1); eqHolm=nan(m,1);
eqRaw(allComputed)=double(globalP(allComputed)<alpha);
eqHolm(allComputed)=double(holmP(allComputed)<alpha);
col1=pn(pairs(:,1)); col2=pn(pairs(:,2));
name1=labels(pairs(:,1)); name2=labels(pairs(:,2));
summary_table=table(col1(:),col2(:),name1(:),name2(:),counts, ...
    repmat(-c,m,1),repmat(c,m,1),allComputed,globalP,eqRaw,holmP,eqHolm, ...
    'VariableNames',{'Parameter1','Parameter2','ParameterName1','ParameterName2', ...
    'N','DzLowerBound','DzUpperBound','AllConditionsComputed','GlobalP', ...
    'EquivalentUnadjusted','GlobalPHolm','EquivalentHolm'});
pm=nan(np); ph=nan(np); em=nan(np);
for h=1:m
    a=pairs(h,1); b=pairs(h,2);
    pm(a,b)=globalP(h); pm(b,a)=globalP(h);
    ph(a,b)=holmP(h); ph(b,a)=holmP(h);
    em(a,b)=eqHolm(h); em(b,a)=eqHolm(h);
end
metadata=struct('Alpha',alpha,'DzBounds',[-c c], ...
    'MeanConfidenceLevel',1-2*alpha,'HolmFamilySize',m, ...
    'ParameterMin',minima,'ParameterMax',maxima, ...
    'Scaling','(X - ParameterMin)/(ParameterMax - ParameterMin)', ...
    'MissingDataRule','Complete cases per parameter pair across four conditions');
condition_table.Properties.UserData=metadata;
summary_table.Properties.UserData=metadata;
details=metadata;
details.ParameterOrder=pn; details.ParameterMax=maxima;
details.ParameterNames=labels; details.ConditionNames=conditions;
details.PairIndices=pairs; details.ValidSubjectsByPair=validSubjects;
details.ConditionP=conditionP; details.GlobalPMatrix=pm;
details.P_Holm_Matrix=ph; details.EquivalentHolmMatrix=em;
fprintf('Parameter levels: %d pairs, four conditions per pair, dz bounds [%g,%g].\n',m,-c,c);
fprintf('GlobalP=max(four TOST p-values); Holm across %d parameter pairs.\n',m);
end

% function ok=positiveScalar(x)
% ok=isnumeric(x)&&isreal(x)&&isscalar(x)&&isfinite(x)&&x>0;
% end


function [pair_table, summary_table] = equivalence_four_conditions_sd(ots, pls, otn, pln, varargin)
%EQUIVALENCE_FOUR_CONDITIONS_SD Paired equivalence of four conditions.
% [pairs, summary] = equivalence_four_conditions_sd(ots,pls,otn,pln,...
%     'SDMultiplier',1,'ParameterOrder',[6 1:5 7 8],'Alpha',0.05)
%
% Each input: subjects x parameters, with identical subject/column ordering.
% For every parameter and condition pair, test -c < mu_D/sigma_D < c.
% Uses noncentral t boundary distributions, not plug-in raw-margin TOST.
% Assumes independent subjects and normally distributed paired differences.
% Estimated raw margins +/-c*SD_D are reported for interpretation only.
% A parameter is globally equivalent only if all six pair tests pass.
% Holm correction across parameter-level global p-values is also reported.
% Pair p-values/flags are unadjusted components of the global test, not
% separately multiplicity-controlled claims. Invalid tests have NaN flags.
% Requires Statistics and Machine Learning Toolbox; no other custom code.

ip = inputParser;
addParameter(ip,'SDMultiplier',1,@(x) isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>0);
addParameter(ip,'Alpha',0.05,@(x) isnumeric(x)&&isscalar(x)&&isfinite(x)&&x>0&&x<0.5);
addParameter(ip,'ParameterOrder',[6 1:5 7 8], ...
    @(x) isnumeric(x)&&isvector(x)&&all(isfinite(x))&&all(x>=1)&&all(x==fix(x)));
parse(ip,varargin{:});
c = ip.Results.SDMultiplier;
alpha = ip.Results.Alpha;
pn = ip.Results.ParameterOrder(:)';
inputs = {ots,pls,otn,pln};
for j = 1:4
    validateattributes(inputs{j},{'numeric'},{'real','2d'});
end
assert(isequal(size(ots),size(pls),size(otn),size(pln)), ...
    'All four matrices must have the same dimensions.');
assert(max(pn)<=size(ots,2) && numel(unique(pn))==numel(pn), ...
    'ParameterOrder must contain unique, existing column indices.');
names = {'OTS','PLS','OTN','PLN'};
comb = nchoosek(1:4,2);
np = numel(pn);
rows = cell(6*np,16);
counts = zeros(np,1);
global_p = nan(np,1);
global_eq = nan(np,1);
complete = false(np,1);

for p = 1:np
    X = double([ots(:,pn(p)),pls(:,pn(p)),otn(:,pn(p)),pln(:,pn(p))]);
    X = X(all(isfinite(X),2),:); % same complete subjects for all six pairs
    n = size(X,1);
    counts(p) = n;
    ps = nan(6,1);
    for k = 1:6
        mu=NaN; sd=NaN; dz=NaN; lo=NaN; hi=NaN;
        pL=NaN; pU=NaN; pt=NaN; eq=NaN;
        status='Too few complete subjects';
        if n >= 2
            d = X(:,comb(k,1))-X(:,comb(k,2));
            mu = mean(d);
            sd = std(d,0);
            status='Zero or nonfinite SD: standardized test unavailable';
            if isfinite(mu) && isfinite(sd) && sd>0
                dz = mu/sd;
                t = dz*sqrt(n);
                % At dz=-c or +c, noncentrality = -c*sqrt(n) or +c*sqrt(n).
                pL = nctcdf(t,n-1,-c*sqrt(n),'upper');
                pU = nctcdf(t,n-1, c*sqrt(n));
                pt = max(pL,pU);
                % Raw CI is descriptive; do not compare it with +/-c*sd
                % to reproduce this standardized-effect decision.
                halfwidth = tinv(1-alpha,n-1)*sd/sqrt(n);
                lo = mu-halfwidth;
                hi = mu+halfwidth;
                if isfinite(pL) && isfinite(pU)
                    eq = double(pt<alpha);
                    status='Computed';
                else
                    pt=NaN;
                    status='Noncentral t calculation unavailable';
                end
            end
        end
        ps(k) = pt;
        rows((p-1)*6+k,:) = {pn(p),names{comb(k,1)},names{comb(k,2)}, ...
            n,mu,sd,dz,-c*sd,c*sd,lo,hi,pL,pU,pt,eq,status};
    end
    complete(p) = all(isfinite(ps));
    if complete(p)
        global_p(p) = max(ps); % intersection-union across six pairs
        global_eq(p) = double(global_p(p)<alpha);
    end
end

pair_table = cell2table(rows,'VariableNames', ...
    {'ParameterColumn','Condition1','Condition2','N','MeanDifference', ...
     'SD_Difference','Dz','EstimatedRawLower','EstimatedRawUpper', ...
     'RawCI_Lower_Unadjusted','RawCI_Upper_Unadjusted', ...
     'P_Lower','P_Upper','P_TOST_Dz','EquivalentUnadjusted','Status'});

% Preserve the whole planned family, including any uncomputable parameters.
pw = global_p;
pw(~isfinite(pw))=1;
[ps,ix] = sort(pw);
adj = min(1,cummax((np:-1:1)'.*ps));
global_holm = nan(np,1);
global_holm(ix)=adj;
global_holm(~complete)=NaN;
eq_holm=nan(np,1);
eq_holm(complete)=double(global_holm(complete)<alpha);
summary_table = table(pn(:),counts,repmat(c,np,1),complete,global_p, ...
    global_eq,global_holm,eq_holm,'VariableNames', ...
    {'ParameterColumn','N','SDMultiplier','AllPairsComputed','GlobalP', ...
     'EquivalentUnadjusted','GlobalPHolm','EquivalentHolm'});
pair_table.Properties.UserData = struct('Alpha',alpha, ...
    'RawConfidenceLevel',1-2*alpha,'DzBounds',[-c c]);
fprintf('Standardized paired equivalence: -%.3g < mu_D/sigma_D < %.3g.\n',c,c);
fprintf('Six pairs per parameter; Holm adjustment across %d parameters.\n',np);
disp(summary_table);
end
