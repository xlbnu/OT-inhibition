
%% cohort1: behavior expriment
[sub_ots,sub_otn]=loadData('E:\xianliang\matlab_m\social_decision_m\data\OT_EXP\OT_data');
[sub_pls,sub_pln]=loadData('E:\xianliang\matlab_m\social_decision_m\data\OT_EXP\PL_data');
%
[ots_bhv,otn_bhv]=getOTdata(sub_ots,sub_otn);
[pls_bhv,pln_bhv]=getOTdata(sub_pls,sub_pln);
[ots_bhv,otn_bhv,pls_bhv,pln_bhv]=correctSubOrder(ots_bhv,otn_bhv,pls_bhv,pln_bhv);

%% cohort2: MEG expriment
[sub_ots,sub_otn]=loadData('E:\xianliang\matlab_m\social_decision_m\data\OT-N2-MEG\behavior\OT_data');
[sub_pls,sub_pln]=loadData('E:\xianliang\matlab_m\social_decision_m\data\OT-N2-MEG\behavior\PL_data');
%
[ots_meg,otn_meg]=getOTdata(sub_ots,sub_otn);
[pls_meg,pln_meg]=getOTdata(sub_pls,sub_pln);
[ots_meg,otn_meg,pls_meg,pln_meg]=correctSubOrder(ots_meg,otn_meg,pls_meg,pln_meg);

%% cohort3: MRS expriment
sub_MRS0=getMRSNumber;
[sub_ots,sub_otn]=loadData('E:\xianliang\matlab_m\social_decision_m\data\OT-N2-MRS\OT_data');
[sub_pls,sub_pln]=loadData('E:\xianliang\matlab_m\social_decision_m\data\OT-N2-MRS\PL_data');
%
[ots_mrs,otn_mrs]=getOTdata(sub_ots,sub_otn,sub_MRS0);
[pls_mrs,pln_mrs]=getOTdata(sub_pls,sub_pln,sub_MRS0);
[ots_mrs,otn_mrs,pls_mrs,pln_mrs]=correctSubOrder(ots_mrs,otn_mrs,pls_mrs,pln_mrs);

%% merge three cohort experiments' data
ots_data=getConjunctionData({ots_bhv;ots_meg;ots_mrs});
pls_data=getConjunctionData({pls_bhv;pls_meg;pls_mrs});
otn_data=getConjunctionData({otn_bhv;otn_meg;otn_mrs});
pln_data=getConjunctionData({pln_bhv;pln_meg;pln_mrs});

save('E:\xianliang\matlab_m\social_decision_m\data\OT_all\BHV_MEG_MRS_data_all\OT_data.mat','ots_data','pls_data','otn_data','pln_data');


%% functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [sub_social, sub_infor] = loadData(ndir)
% Load each participant folder Social / nonSocial data.

folders = dir(ndir);
folders = folders([folders.isdir] & ~ismember({folders.name}, {'.', '..'}));
sub_social = cell(numel(folders), 1);
sub_infor = cell(numel(folders), 1);

for i = 1:numel(folders)
    folderName = folders(i).name;
    folderPath = fullfile(folders(i).folder, folderName);
    subjectID = regexprep(folderName, '-\d+$', '');
    sessionID = regexprep(folderName, '-(?=\d+$)', '');
    files = dir(folderPath);
    files = files(~[files.isdir]);
    socialData = {};
    nonsocialData = {};

    for j = 1:numel(files)
        [~, filename] = fileparts(files(j).name);
        if ~contains(filename, 'Social')
            continue;
        end

        % Load behavioral data, and record the experiment time from the filename.
        loaded = load(fullfile(folderPath, files(j).name), 'savedata');
        data = loaded.savedata.data;
        time = regexp(filename, '\d{2};\d{2};\d{2}', 'match', 'once');
        data.session_time = [extractBefore(filename, 6), ';', time];

        % Check participant names, Attach the experiment-order number from the folder.
        isNonsocial = contains(filename, 'nonSocial');
        if isNonsocial
            name = data.subsName;
        else
            name = data.subsName{1};
        end
        if ~isequal(name, subjectID)
            error('subject name is not equal');
        end

        % Preserve the original output format: social participant name is cell, nonsocial is a character array.
        if isNonsocial
            data.subsName = sessionID;
            nonsocialData{end + 1} = data;
        else
            data.subsName{1} = sessionID;
            socialData{end + 1} = data;
        end
    end

    sub_social{i} = struct('expdata', socialData);
    sub_infor{i} = struct('expdata', nonsocialData);
end
end




function [ots_data, otn_data] = getOTdata(sub_ots, sub_otn, sub_MRS0)
% Merge experimental data for each participant, Compute choice, confidence and log-odds measures.

ots_data = struct([]);
otn_data = struct([]);

for i = 1:numel(sub_ots)
    % Social: Merge data, Reconstruct confidence using key responses and left/right mappings.
    social = mergeSessions(sub_ots{i}, true);
    for phase = 1:2
        keyField = sprintf('ckey%d', phase);
        confField = sprintf('conf%d', phase);
        keys = social.(keyField);
        keys(keys < 1 | keys > 8) = NaN;
        social.(keyField) = keys;
        keys(social.rlo == 1) = 9 - keys(social.rlo == 1);
        social.(confField) = keys;
    end

    social.conflict = social.answer1 ~= social.o1answer1;
    otherProbability = social.o1conf1 / 16 + 0.5 - 1/32;
    social = computeDecisionVariables(social, otherProbability);

    % Extract participant name, experiment order and MRS identifier.
    [name, order] = subjectInfo(social.subsName{1});
    social.subsName = {name};
    social.expOrder = order;
    social.mrs_num = sub_MRS0{ismember(sub_MRS0(:, 2), name), 3};
    ots_data = [ots_data; social];

    % nonSocial: Use other-person evidence directly from inforate.
    nonsocial = mergeSessions(sub_otn{i}, false);
    nonsocial.conflict = nonsocial.answer1 ~= nonsocial.infoans;
    nonsocial = computeDecisionVariables(nonsocial, nonsocial.inforate);

    % Missing nonSocial participant name:, use the corresponding Social identity and order.
    if isfield(nonsocial, 'subsName')
        sessions = [sub_otn{i}.expdata];
        [nonName, nonOrder] = subjectInfo(sessions(1).subsName);
        if numel(nonName) < 2
            nonName = name;
        end
    else
        nonName = name;
        nonOrder = order;
    end
    nonsocial.subsName = nonName;
    nonsocial.expOrder = nonOrder;
    nonsocial.mrs_num = sub_MRS0{ismember(sub_MRS0(:, 2), nonName), 3};
    otn_data = [otn_data; nonsocial];
end
end

function data = mergeSessions(subject, isSocial)
% Merge sessions in their original order; Preserve the original transpose rules for wide arrays.
sessions = [subject.expdata];
if isSocial && isfield(sessions, 'dsOrder')
    sessions = rmfield(sessions, 'dsOrder');
end

names = fieldnames(sessions);
data = struct();
for f = 1:numel(names)
    field = names{f};
    for j = 1:numel(sessions)
        value = sessions(j).(field);
        if size(value, 2) > 4 && (isSocial || ~ischar(value))
            sessions(j).(field) = value';
        end
    end
    data.(field) = cat(1, sessions.(field));
end
end

function data = computeDecisionVariables(data, otherProbability)
% On congruent trials, if confidence did not decrease after switching, restore the initial choice using the original rule.
restoreChoice = ~data.conflict & data.answer1 ~= data.answer2 & data.conf1 <= data.conf2;
data.answer2(restoreChoice) = data.answer1(restoreChoice);
data.isc2 = data.answer2 == data.dir;
data.ischange = data.answer1 ~= data.answer2;

% Convert own confidence and other-person evidence to log-odds.
logOdds = @(p) log(p ./ (1 - p));
data.lo_s1 = logOdds(data.conf1 / 16 + 0.5 - 1/32);
data.lo_s2 = logOdds(data.conf2 / 16 + 0.5 - 1/32);
data.lo_o1 = logOdds(otherProbability);
data.lo_plus = data.lo_s1 + data.lo_o1;
data.lo_minus = data.lo_s1 - data.lo_o1;

% When staying, own evidence is chosen; When switching, other-person evidence is chosen.
data.lo_chosen = data.lo_s1;
data.lo_unchosen = data.lo_o1;
changed = data.ischange == 1;
data.lo_chosen(changed) = data.lo_o1(changed);
data.lo_unchosen(changed) = data.lo_s1(changed);
data.lo_cplus = data.lo_chosen + data.lo_unchosen;
data.lo_cminus = data.lo_chosen - data.lo_unchosen;
end

function [name, order] = subjectInfo(sessionName)
% Digits in the name specify experiment order; Default to the first session when no digits are present.
if iscell(sessionName)
    sessionName = sessionName{1};
end
digits = regexp(sessionName, '\d+', 'match', 'once');
order = str2double(digits);
if isempty(digits)
    order = 1;
end
name = lower(regexprep(sessionName, '\d', ''));
end



% make sure the subjects' order same between ot-data and pl-data
function [ots_data,otn_data,pls_data,pln_data]=correctSubOrder(ots_data,otn_data,pls_data,pln_data)
ot_name=cell(length(ots_data),1);ot1_name=cell(length(otn_data),1);
pl_name=cell(length(pls_data),1);pl1_name=cell(length(pln_data),1);
for i=1:length(ots_data)
ot_name{i,1}=ots_data(i).subsName{1};ot1_name{i,1}=otn_data(i).subsName;
pl_name{i,1}=pls_data(i).subsName{1};pl1_name{i,1}=pln_data(i).subsName;
end
if ~isequal(ot_name,ot1_name) | ~isequal(pl_name,pl1_name)
error('social & nonsocial: subjects name is not in the correct order');
end
if length(ot_name)~=length(pl_name)
error('ot-data & pl-data: the number between ot-data and pl-data are not same');
end

[~,idx]=ismember(ot_name,pl_name);
pls_data=pls_data(idx,:);
pln_data=pln_data(idx,:);
pl_name1=cell(length(pls_data),1);pl1_name1=cell(length(pln_data),1);
for i=1:length(pls_data)
pl_name1{i,1}=pls_data(i).subsName{1};pl1_name1{i,1}=pln_data(i).subsName;
end
if ~isequal(ot_name,pl_name1) | ~isequal(ot1_name,pl1_name1) | ~isequal(pl_name1,pl1_name1)
error('something wrong');
end
%%
% Sort the two strings in each row(Standardize formatting)
C=cat(1,ots_data.subsName);sortedStr=cell(size(C,1),size(C,2));
for i = 1:size(C, 1)
    rowStrings = C(i, :);     % Extract the two strings in the current row
    sortedRow = sort(rowStrings);  % Sort lexicographically
    sortedStr(i, :) = sortedRow;   % Store the sorted result
end
% Create a unique key for the sorted strings(Avoid comparing directly using'rows'parameters cell)
[~, ~, groupID] = unique(...
    strcat(sortedStr(:, 1), '_', sortedStr(:, 2)));  % Join strings using a delimiter
for i=1:length(ots_data)
ots_data(i).pair_num=groupID(i);otn_data(i).pair_num=groupID(i);
pls_data(i).pair_num=groupID(i);pln_data(i).pair_num=groupID(i);
end

end


function sub_MRS0=getMRSNumber
path_date_data='E:\xianliang\matlab_m\social_decision_m\data\OT-N2-MRS\sendout\dataBydate';
n_data=dir(path_date_data);
sub_MRS=[];
for i=1:length(n_data)-2
    load([path_date_data,'\',n_data(i+2).name]);
    sub_MRS=cat(1,sub_MRS,sub_data(:,2:3));    
end

[~,idx]=unique(sub_MRS(:,1));
sub_MRS0=sub_MRS(idx,1:2);
sub_MRS0([1],:)=[];
for i=1:length(sub_MRS0)
sub_MRS0{i,3}=str2double(regexprep(sub_MRS0{i,1},'S0',''));
if sub_MRS0{i,3}>9
sub_MRS0{i,4}=['S',mat2str(sub_MRS0{i,3})];
else
sub_MRS0{i,4}=['S0',mat2str(sub_MRS0{i,3})];    
end
end
end


function ots_all = getConjunctionData(dataCell)
    % Example input: { {structArray1, structArray2}, {structArray3} }
    if ~iscell(dataCell)
        error('输入必须为 cell 数组。');
    end

    % Recursively extract all nonempty structure arrays
    structList = collectStructs(dataCell);

    if isempty(structList)
        ots_all = struct([]);
        return;
    end

    % Get fields shared by all structures
    commonFields = fieldnames(structList{1});
    for i = 2:numel(structList)
        commonFields = intersect(commonFields, ...
            fieldnames(structList{i}), 'stable');
    end

    if isempty(commonFields)
        error('各结构体之间没有共同字段，无法合并。');
    end

    % Retain shared fields, Standardize field order and array orientation
    for i = 1:numel(structList)
        s = structList{i};
        s = rmfield(s, setdiff(fieldnames(s), commonFields));
        s = orderfields(s, commonFields);
        structList{i} = s(:);
    end

    ots_all = vertcat(structList{:});
end

function structList = collectStructs(data)
    structList = {};

    if iscell(data)
        for i = 1:numel(data)
            subList = collectStructs(data{i});
            structList = [structList; subList]; %#ok<AGROW>
        end
    elseif isstruct(data)
        if ~isempty(data)
            structList = {data};
        end
    elseif ~isempty(data)
        error('cell 中只能包含结构体、子 cell 或空值。');
    end
end