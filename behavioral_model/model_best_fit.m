
% Run this script once per selected analysis; fitting must finish first.
codeRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(codeRoot,'utilities'));
if ~exist('bestFitAnalysis','var') || isempty(bestFitAnalysis)
    bestFitAnalysis = 'full'; % 'full', 'two_lambda', or 'comparison'
end
assert(ismember(bestFitAnalysis,{'full','two_lambda','comparison'}),'Invalid bestFitAnalysis.');
if exist('bestFitInputDir','var') && ~isempty(bestFitInputDir)
    raw_dir = bestFitInputDir;
elseif strcmp(bestFitAnalysis,'full')
    raw_dir = latest_run_local(fullfile(codeRoot,'analysis_outputs','model_fitting','full'));
elseif strcmp(bestFitAnalysis,'two_lambda')
    raw_dir = latest_run_local(fullfile(codeRoot,'analysis_outputs','model_fitting','two_lambda'));
else
    raw_dir = latest_run_local(fullfile(codeRoot,'analysis_outputs','model_comparison'));
end
% Set bestFitInputDir before running to select a specific RUN_TAG directory.
outputDir = fullfile(codeRoot,'analysis_outputs','model_results');
if ~isfolder(outputDir), mkdir(outputDir); end

if strcmp(bestFitAnalysis,'full')
[os_fit1,ps_fit1,on_fit1,pn_fit1] = read_repeats_local(raw_dir);
%% FOR FULL MODEL (NINE PARAMETER)
os_fit2=cell(length(os_fit1),1);ps_fit2=cell(length(ps_fit1),1);
on_fit2=cell(length(on_fit1),1);pn_fit2=cell(length(pn_fit1),1);
for i=1:length(os_fit1)
    os_valid0=getNoExtrmeFitParams(os_fit1{i});
    os_fit2{i,1}=sortrows(os_fit1{i}(os_valid0,:),10);
    ps_valid0=getNoExtrmeFitParams(ps_fit1{i});
    ps_fit2{i,1}=sortrows(ps_fit1{i}(ps_valid0,:),10);
    on_valid0=getNoExtrmeFitParams(on_fit1{i});
    on_fit2{i,1}=sortrows(on_fit1{i}(on_valid0,:),10);
    pn_valid0=getNoExtrmeFitParams(pn_fit1{i});
    pn_fit2{i,1}=sortrows(pn_fit1{i}(pn_valid0,:),10);
end




%
ots_minb1=[];pls_minb1=[];otn_minb1=[];pln_minb1=[];
for i=1:length(os_fit2)

    ots_tmp.optParams=os_fit2{i}(1,1:9);
    ots_tmp.bic=os_fit2{i}(1,10);
    ots_tmp.nll=os_fit2{i}(1,11);
    ots_minb1=cat(1,ots_minb1,ots_tmp);

    pls_tmp.optParams=ps_fit2{i}(1,1:9);
    pls_tmp.bic=ps_fit2{i}(1,10);
    pls_tmp.nll=ps_fit2{i}(1,11);
    pls_minb1=cat(1,pls_minb1,pls_tmp);   

    otn_tmp.optParams=on_fit2{i}(1,1:9);
    otn_tmp.bic=on_fit2{i}(1,10);
    otn_tmp.nll=on_fit2{i}(1,11);
    otn_minb1=cat(1,otn_minb1,otn_tmp);

    pln_tmp.optParams=pn_fit2{i}(1,1:9);
    pln_tmp.bic=pn_fit2{i}(1,10);
    pln_tmp.nll=pn_fit2{i}(1,11);
    pln_minb1=cat(1,pln_minb1,pln_tmp);    
    
end

%
ots_xx1=cat(1,ots_minb1.optParams);
pls_xx1=cat(1,pls_minb1.optParams);
otn_xx1=cat(1,otn_minb1.optParams);
pln_xx1=cat(1,pln_minb1.optParams);
%
pn=[7 1:6 8 9];p_params=nan(9,1);
for i=1:9
[~,p_params(i,1)]=ttest(ots_xx1(:,pn(i))-otn_xx1(:,pn(i))-(pls_xx1(:,pn(i))-pln_xx1(:,pn(i))));
end

%
x_pos = [0.5 1.7 3 4.2];
figureX1([28 10]);
p_ylim=[-0 0 0 0 0 0 0 0 0;1.5 1.5 1.5 15 15 25 0.2 1500 1500]';
paramNames = {'v_0', 'b_1/b_0', 'b_2/b_0', 'k1', 'k2', 'σ', 'λ', '\tau', 't_0'};
pn=[7 1:6 8 9];
tiledlayout(3,3);
for i=1:9
    nexttile;
    hold on;
    figureTmp1([ots_xx1(:,pn(i)),pls_xx1(:,pn(i)),otn_xx1(:,pn(i)),pln_xx1(:,pn(i))]);   
    xx1=[ots_xx1(:,pn(i)),pls_xx1(:,pn(i)),otn_xx1(:,pn(i)),pln_xx1(:,pn(i))];
% [~,p]=ttest(xx1(:,1)-xx1(:,3)-(xx1(:,2)-xx1(:,4)))
    hold off;box off;ylabel(paramNames{pn(i)},'Interpreter','tex');set(gca,'XTick',x_pos,'XTickLabel',{'OT','PL','OT','PL'});xlim([-0.2 4.9]); gcaf1(gca);
end


save(fullfile(outputDir,'model_fit_parameter_85.mat'),'ots_minb1','pls_minb1','otn_minb1','pln_minb1');
end

if strcmp(bestFitAnalysis,'comparison')
%% model compare： bic \ aic \ nll

rep_dir=dir(fullfile(raw_dir,'ots_data','*_*'));
bic_os=nan(length(rep_dir),count_repeats_local(raw_dir),85);nll_os=bic_os;aic_os=bic_os;
for i=1:length(rep_dir)
    rep_dir1=dir(fullfile(rep_dir(i).folder,rep_dir(i).name,'repeat_*'));
    assert(~isempty(rep_dir1),'No repeat folders for comparison model.');
    for j=1:length(rep_dir1)
        assert(isfile(fullfile(rep_dir1(j).folder,rep_dir1(j).name,'results_all.mat')), ...
            'Incomplete comparison run: results_all.mat is missing.');
        if isfile(fullfile(rep_dir1(j).folder,rep_dir1(j).name,'results_all.mat'))
            xx_tmp=load(fullfile(rep_dir1(j).folder,rep_dir1(j).name,'results_all.mat'),'results');
            os_fit{i,j}=[cat(1,xx_tmp.results.optParams_FullArray),cat(1,xx_tmp.results.bic),cat(1,xx_tmp.results.nll)];
            bic_os(i,j,:)=cat(1,xx_tmp.results.bic);
            nll_os(i,j,:)=cat(1,xx_tmp.results.nll);
            aic_os(i,j,:)=cat(1,xx_tmp.results.nll).*2+2*(10-i);            
        end
    end
    assert(exist('xx_tmp','var') == 1,'No completed fits for comparison model.');
    os_fit{i,j+1}=xx_tmp.results(1).modelName;
end

rep_dir=dir(fullfile(raw_dir,'pls_data','*_*'));
bic_ps=nan(length(rep_dir),count_repeats_local(raw_dir),85);nll_ps=bic_ps;aic_ps=bic_ps;
for i=1:length(rep_dir)
    rep_dir1=dir(fullfile(rep_dir(i).folder,rep_dir(i).name,'repeat_*'));
    assert(~isempty(rep_dir1),'No repeat folders for comparison model.');
    for j=1:length(rep_dir1)
        assert(isfile(fullfile(rep_dir1(j).folder,rep_dir1(j).name,'results_all.mat')), ...
            'Incomplete comparison run: results_all.mat is missing.');
        if isfile(fullfile(rep_dir1(j).folder,rep_dir1(j).name,'results_all.mat'))
            xx_tmp=load(fullfile(rep_dir1(j).folder,rep_dir1(j).name,'results_all.mat'),'results');
            ps_fit{i,j}=[cat(1,xx_tmp.results.optParams_FullArray),cat(1,xx_tmp.results.bic),cat(1,xx_tmp.results.nll)];
            bic_ps(i,j,:)=cat(1,xx_tmp.results.bic);
            nll_ps(i,j,:)=cat(1,xx_tmp.results.nll);
            aic_ps(i,j,:)=cat(1,xx_tmp.results.nll).*2+2*(10-i);
        end
    end
    assert(exist('xx_tmp','var') == 1,'No completed fits for comparison model.');
    ps_fit{i,j+1}=xx_tmp.results(1).modelName;
end
%
rep_dir=dir(fullfile(raw_dir,'otn_data','*_*'));
bic_on=nan(length(rep_dir),count_repeats_local(raw_dir),85);nll_on=bic_on;aic_on=bic_on;
for i=1:length(rep_dir)
    rep_dir1=dir(fullfile(rep_dir(i).folder,rep_dir(i).name,'repeat_*'));
    assert(~isempty(rep_dir1),'No repeat folders for comparison model.');
    for j=1:length(rep_dir1)
        assert(isfile(fullfile(rep_dir1(j).folder,rep_dir1(j).name,'results_all.mat')), ...
            'Incomplete comparison run: results_all.mat is missing.');
        if isfile(fullfile(rep_dir1(j).folder,rep_dir1(j).name,'results_all.mat'))
            xx_tmp=load(fullfile(rep_dir1(j).folder,rep_dir1(j).name,'results_all.mat'),'results');
            on_fit{i,j}=[cat(1,xx_tmp.results.optParams_FullArray),cat(1,xx_tmp.results.bic),cat(1,xx_tmp.results.nll)];
            bic_on(i,j,:)=cat(1,xx_tmp.results.bic);
            nll_on(i,j,:)=cat(1,xx_tmp.results.nll);
            aic_on(i,j,:)=cat(1,xx_tmp.results.nll).*2+2*(10-i);
        end
    end
    assert(exist('xx_tmp','var') == 1,'No completed fits for comparison model.');
    on_fit{i,j+1}=xx_tmp.results(1).modelName;
end
%
rep_dir=dir(fullfile(raw_dir,'pln_data','*_*'));
bic_pn=nan(length(rep_dir),count_repeats_local(raw_dir),85);nll_pn=bic_pn;aic_pn=bic_pn;
for i=1:length(rep_dir)
    rep_dir1=dir(fullfile(rep_dir(i).folder,rep_dir(i).name,'repeat_*'));
    assert(~isempty(rep_dir1),'No repeat folders for comparison model.');
    for j=1:length(rep_dir1)
        assert(isfile(fullfile(rep_dir1(j).folder,rep_dir1(j).name,'results_all.mat')), ...
            'Incomplete comparison run: results_all.mat is missing.');
        if isfile(fullfile(rep_dir1(j).folder,rep_dir1(j).name,'results_all.mat'))
            xx_tmp=load(fullfile(rep_dir1(j).folder,rep_dir1(j).name,'results_all.mat'),'results');
            pn_fit{i,j}=[cat(1,xx_tmp.results.optParams_FullArray),cat(1,xx_tmp.results.bic),cat(1,xx_tmp.results.nll)];
            bic_pn(i,j,:)=cat(1,xx_tmp.results.bic);
            nll_pn(i,j,:)=cat(1,xx_tmp.results.nll);
            aic_pn(i,j,:)=cat(1,xx_tmp.results.nll).*2+2*(10-i);
        end
    end
    assert(exist('xx_tmp','var') == 1,'No completed fits for comparison model.');
    pn_fit{i,j+1}=xx_tmp.results(1).modelName;
end

%
mc_bic.os=squeeze(min(bic_os(:,:,:),[],2))';
mc_bic.ps=squeeze(min(bic_ps(:,:,:),[],2))';
mc_bic.on=squeeze(min(bic_on(:,:,:),[],2))';
mc_bic.pn=squeeze(min(bic_pn(:,:,:),[],2))';

mc_nll.os=squeeze(min(nll_os(:,:,:),[],2))';
mc_nll.ps=squeeze(min(nll_ps(:,:,:),[],2))';
mc_nll.on=squeeze(min(nll_on(:,:,:),[],2))';
mc_nll.pn=squeeze(min(nll_pn(:,:,:),[],2))';

mc_aic.os=squeeze(min(aic_os(:,:,:),[],2))';
mc_aic.ps=squeeze(min(aic_ps(:,:,:),[],2))';
mc_aic.on=squeeze(min(aic_on(:,:,:),[],2))';
mc_aic.pn=squeeze(min(aic_pn(:,:,:),[],2))';

[~,tmp]=min(bic_os(:,:,:),[],2);idx_bic.os=squeeze(tmp);
[~,tmp]=min(bic_ps(:,:,:),[],2);idx_bic.ps=squeeze(tmp);
[~,tmp]=min(bic_on(:,:,:),[],2);idx_bic.on=squeeze(tmp);
[~,tmp]=min(bic_pn(:,:,:),[],2);idx_bic.pn=squeeze(tmp);
%
save(fullfile(outputDir,'model_compare_data.mat'),'mc_bic','mc_aic','mc_nll');















end

if strcmp(bestFitAnalysis,'two_lambda')
[os_fit1,ps_fit1,on_fit1,pn_fit1] = read_repeats_local(raw_dir);
%% FOR bidirectional inhibition MODEL
os_fit2=cell(length(os_fit1),1);ps_fit2=cell(length(ps_fit1),1);
on_fit2=cell(length(on_fit1),1);pn_fit2=cell(length(pn_fit1),1);
for i=1:length(os_fit1)
    os_valid0=getNoExtrmeFitParams2(os_fit1{i});
    os_fit2{i,1}=sortrows(os_fit1{i}(os_valid0,:),11);
    ps_valid0=getNoExtrmeFitParams2(ps_fit1{i});
    ps_fit2{i,1}=sortrows(ps_fit1{i}(ps_valid0,:),11);
    on_valid0=getNoExtrmeFitParams2(on_fit1{i});
    on_fit2{i,1}=sortrows(on_fit1{i}(on_valid0,:),11);
    pn_valid0=getNoExtrmeFitParams2(pn_fit1{i});
    pn_fit2{i,1}=sortrows(pn_fit1{i}(pn_valid0,:),11);
end




%%
ots_minb1=[];pls_minb1=[];otn_minb1=[];pln_minb1=[];
for i=1:length(os_fit2)

    ots_tmp.optParams=os_fit2{i}(1,1:10);
    ots_tmp.bic=os_fit2{i}(1,11);
    ots_tmp.nll=os_fit2{i}(1,12);
    ots_minb1=cat(1,ots_minb1,ots_tmp);

    pls_tmp.optParams=ps_fit2{i}(1,1:10);
    pls_tmp.bic=ps_fit2{i}(1,11);
    pls_tmp.nll=ps_fit2{i}(1,12);
    pls_minb1=cat(1,pls_minb1,pls_tmp);   

    otn_tmp.optParams=on_fit2{i}(1,1:10);
    otn_tmp.bic=on_fit2{i}(1,11);
    otn_tmp.nll=on_fit2{i}(1,12);
    otn_minb1=cat(1,otn_minb1,otn_tmp);

    pln_tmp.optParams=pn_fit2{i}(1,1:10);
    pln_tmp.bic=pn_fit2{i}(1,11);
    pln_tmp.nll=pn_fit2{i}(1,12);
    pln_minb1=cat(1,pln_minb1,pln_tmp);    
    
end

%
ots_xx1=cat(1,ots_minb1.optParams);
pls_xx1=cat(1,pls_minb1.optParams);
otn_xx1=cat(1,otn_minb1.optParams);
pln_xx1=cat(1,pln_minb1.optParams);
%
pn=[7 1:6 8 9];p_params=nan(9,1);
for i=1:9
[~,p_params(i,1)]=ttest(ots_xx1(:,pn(i))-otn_xx1(:,pn(i))-(pls_xx1(:,pn(i))-pln_xx1(:,pn(i))));
end

%
x_pos = [0.5 1.7 3 4.2];
figureX1([28 10]);
p_ylim=[-0 0 0 0 0 0 0 0 0;1.5 1.5 1.5 15 15 25 0.2 1500 1500]';
paramNames = {'v_0', 'b_1/b_0', 'b_2/b_0', 'k1', 'k2', 'σ', 'λ1','λ2', '\tau', 't_0'};
pn=[7 8 1:6 9 10];
tiledlayout(2,5);
for i=1:10
    nexttile;
    hold on;
    figureTmp1([ots_xx1(:,pn(i)),pls_xx1(:,pn(i)),otn_xx1(:,pn(i)),pln_xx1(:,pn(i))]);   
    xx1=[ots_xx1(:,pn(i)),pls_xx1(:,pn(i)),otn_xx1(:,pn(i)),pln_xx1(:,pn(i))];
% [~,p]=ttest(xx1(:,1)-xx1(:,3)-(xx1(:,2)-xx1(:,4)))
    hold off;box off;ylabel(paramNames{pn(i)},'Interpreter','tex');set(gca,'XTick',x_pos,'XTickLabel',{'OT','PL','OT','PL'});xlim([-0.2 4.9]); gcaf1(gca);
end











save(fullfile(outputDir,'twolambda_model_fit_parameter_85.mat'),'ots_minb1','pls_minb1','otn_minb1','pln_minb1');
end

function os_valid0=getNoExtrmeFitParams(os_fit1)
lb = [0.001, 0, 0, 0.001, 0.001, 0.1, 1e-5, 10, 10];
ub = [0.99, 1, 1, 10, 10, 20, 0.15, 1000, 1000]; % Reduceτ,t0 upper bound
os_valid=zeros(size(os_fit1,1),length(lb));
for j=[1:7 8 9]
    os_valid(:,j)=os_fit1(:,j)==lb(j) | os_fit1(:,j)==ub(j);
end
os_valid0=sum(os_valid,2)<=1;
end

function os_valid0=getNoExtrmeFitParams2(os_fit1)
lb = [0.001, 0, 0, 0.001, 0.001, 0.1, 1e-5,1e-5, 10, 10];
ub = [0.99, 1, 1, 10, 10, 20, 0.15,0.15, 1000, 1000]; % Reduceτ,t0 upper bound
os_valid=zeros(size(os_fit1,1),length(lb));
for j=[1:7 8 9 10]
    os_valid(:,j)=os_fit1(:,j)==lb(j) | os_fit1(:,j)==ub(j);
end
os_valid0=sum(os_valid,2)<=1;
end

%%  two pair data,distribution+boxplot+scatter+lineconnect
function figureTmp1(xx1, ot_color)
% figureTmp1 - Plot paired data distributions and add significance markers
% Input parameters: 
%   xx1: Data matrix, Each column is a dataset
%   ot_color: Optional color matrix, Default is 4 colors

% Set default colors
if nargin < 2
    ot_color = [197 90 17; 248 203 173; 84 130 53; 197 224 180] ./ 255;
end

hold on;

% Define graphical element positions
x_pos = [0.5 1.7 3 4.2];

% Compute y Maximum axis value and range
ymax_data = max(xx1, [], 'all');
ymin_data = min(xx1, [], 'all');
y_range = ymax_data - ymin_data;

% Plot kernel-density estimates
for i = 1:2
    % First dataset
    x_tmp = xx1(:, 2*i-1);
    [f, xi] = ksdensity(x_tmp, 'NumPoints', 500);
    f = f ./ max(f) .* 0.5;
    ql0 = max(x_tmp) .* 1.01;
    qu0 = min(x_tmp) .* 0.98;
    xu = find(xi > ql0, 1);
    xl = find(xi < qu0, 1, 'last');
    f(xi > ql0 | xi < qu0) = 0;
    xi(xu) = xi(xu-1);
    xi(xl) = xi(xl+1);
    fill(-f + x_pos(2*i-1), xi, ot_color(2*i-1, :), 'EdgeColor', [1 1 1]);

    % Second dataset
    x_tmp = xx1(:, 2*i);
    [f, xi] = ksdensity(x_tmp, 'NumPoints', 500);
    f = f ./ max(f) .* 0.5;
    ql0 = max(x_tmp) .* 1.01;
    qu0 = min(x_tmp) .* 0.98;
    xu = find(xi > ql0, 1);
    xl = find(xi < qu0, 1, 'last');
    f(xi > ql0 | xi < qu0) = 0;
    xi(xu) = xi(xu-1);
    xi(xl) = xi(xl+1);
    fill(f + x_pos(2*i), xi, ot_color(2*i, :), 'EdgeColor', [1 1 1], 'FaceAlpha', 1);
end

% Plot paired connecting lines
x_pos1 = rand(size(xx1, 1), size(xx1, 2)) * 0.2 + x_pos + [0.15 -0.35 0.15 -0.35] + [0.1 -0.1 0.1 -0.1];
pl1 = plot(x_pos1(:, 1:2)', xx1(:, [1 2])', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.4);
pl2 = plot(x_pos1(:, 3:4)', xx1(:, [3 4])', 'Color', [0.7 0.7 0.7], 'LineWidth', 0.4);
for i = 1:length(pl1)
    pl1(i).Color(4) = 0.5;
    pl2(i).Color(4) = 0.5;
end

% Plot scatter points
sc = scatter(x_pos1, xx1(:, [1 2 3 4]), 'MarkerEdgeColor', [1 1 1], ...
             'MarkerFaceColor', 'flat', 'SizeData', 10, 'MarkerFaceAlpha', 0.7);
for i = 1:4
    sc(i).MarkerFaceColor = ot_color(i, :);
end

% Plot boxplots
boxplot(xx1, 'Positions', x_pos + [0.1 -0.1 0.1 -0.1], ...
        'Color', [1 1 1], 'Widths', 0.15, 'Symbol', '');
boxObj = findobj(gca, 'Tag', 'Box'); % Find box objects
for i = 1:length(boxObj)
    patch(get(boxObj(i), 'XData'), get(boxObj(i), 'YData'), ...
          ot_color(5-i, :), 'FaceAlpha', 0.9);
end
boxplot(xx1, 'Positions', x_pos + [0.1 -0.1 0.1 -0.1], ...
        'color', [0 0 0], 'Widths', 0.15, 'Symbol', '');
RemoveBoxplotlines(gca);

% Compute significance
sig_texts = cell(1, 2);
sig_levels = zeros(1, 2);
for i = 1:2
    pair_idx = 2*i-1:2*i;
    [~, p] = ttest(xx1(:, pair_idx(1)), xx1(:, pair_idx(2)));
    if p < 0.001
        sig_texts{i} = '***';
    elseif p < 0.01
        sig_texts{i} = '**';
    elseif p < 0.05
        sig_texts{i} = '*';
    else
        sig_texts{i} = 'n.s.';
    end
    sig_levels(i) = p;
end

% Between-group comparison(Repeated measures ANOVA)
ranovatbl = p2Ride2Ranova(xx1);
if ranovatbl.pValue(7) < 0.001
    between_sig = '***';
elseif ranovatbl.pValue(7) < 0.01
    between_sig = '**';
elseif ranovatbl.pValue(7) < 0.05
    between_sig = '*';
else
    between_sig = 'n.s.';
end

% Set y axis range
% First set an initial range
y_min_plot = ymin_data - 0.1 * y_range;
y_max_plot = ymax_data + 0.3 * y_range;  % Leave space for significance markers
ylim([y_min_plot, y_max_plot]);

% Custom y axis ticks, and return the adjusted y axis range
[new_y_min, new_y_max, increment] = setCustomYTicks(gca, y_min_plot, y_max_plot, xx1);

% Recalculate significance-marker positions(Based on the adjusted y axis range)
% recalculate y axis range
new_y_range = new_y_max - new_y_min;

% Compute vertical positions of significance markers
% First level: Within-group comparison(Paired comparison)
sig_level1_y = new_y_max - 0.20 * new_y_range;
% Second level: Between-group comparison
sig_level2_y = new_y_max - 0.10 * new_y_range;
% Text position for paired-data significance markers(Between the two horizontal lines)
sig_level_mid_y = (sig_level1_y + sig_level2_y) / 2;
% Text height
sig_text_offset = 0.01 * new_y_range;

% Add within-group significance markers
for i = 1:2
    x_pair = x_pos(2*i-1:2*i);
    
    % Draw horizontal lines, Set line width to 0.5(Default width)
    plot(x_pair, [sig_level1_y, sig_level1_y], 'k-', 'LineWidth', 0.5);
    
    % Add significance text
    x_c = mean(x_pair);
    text_x = x_c;
    text_y = sig_level_mid_y;
    
    text(text_x, text_y, sig_texts{i}, ...
         'HorizontalAlignment', 'center', ...
         'VerticalAlignment', 'middle', ...  % Use middle Alignment, Center the text between the two horizontal lines
         'FontSize', 8, ...
         'FontName', 'Arial');
end

% Add between-group significance markers
x_group1 = mean(x_pos(1:2));
x_group2 = mean(x_pos(3:4));
x_line = [x_group1, x_group2];

% Draw between-group comparison lines, Set line width to 0.5(Default width)
plot(x_line, [sig_level2_y, sig_level2_y], 'k-', 'LineWidth', 0.5);

% Add between-group significance text
x_c = mean([x_group1, x_group2]);
text_x = x_c;
text_y = sig_level2_y + sig_text_offset;

% Between-group significance markers use bottom Alignment, Keep parallel to the bottom of within-group markers
% First create temporary text to measure text height
temp_text = text(0, 0, 'n.s.', 'FontSize', 8, 'FontName', 'Arial', 'Visible', 'off');
text_extent = get(temp_text, 'Extent');
text_height = text_extent(4);
delete(temp_text);

% Compute between-group marker text position, Align its bottom with the within-group marker bottom
% Within-group markers use middle Alignment, Bottom at sig_level_mid_y - text_height/2
% Between-group markers use bottom Alignment, y coordinate should match the within-group marker bottom
group_bottom_y = sig_level_mid_y - text_height/2;
% To keep between-group markers above the line, use the original position
% Keep the original position unchanged here
text(text_x, text_y, between_sig, ...
     'HorizontalAlignment', 'center', ...
     'VerticalAlignment', 'bottom', ...  % Use bottom Alignment
     'FontSize', 8, ...
     'FontName', 'Arial');

% Improve plot appearance
box off;
gcaf1(gca);
set(gca, 'TickDir', 'out');
grid off;

% Set x axis labels
set(gca, 'XTick', [mean(x_pos(1:2)), mean(x_pos(3:4))], ...
         'XTickLabel', {'Social', 'Nonsocial'});

end

function [new_y_min, new_y_max, increment] = setCustomYTicks(ax, y_min_plot, y_max_plot, data)
% Set custom y axis ticks
% Input parameters: 
%   ax: axes handle
%   y_min_plot: figure y Axis minimum
%   y_max_plot: figure y Axis maximum
%   data: Original data, Used to determine whether to include 0 tick
% Output parameters: 
%   new_y_min: Adjusted y Axis minimum
%   new_y_max: Adjusted y Axis maximum
%   increment: Tick interval used

% Compute data minimum and maximum
ymin_data = min(data, [], 'all');
ymax_data = max(data, [], 'all');

% Determine the target number of ticks
desired_ticks = 5;

% Compute the optimal interval
range_total = y_max_plot - y_min_plot;
if range_total <= 0
    range_total = 1;  % Prevent division by zero
end
ideal_increment = range_total / desired_ticks;

% Improved interval-selection logic, Adapt to different data ranges
if ideal_increment <= 0
    increment = 0.1;
else
    % Compute logarithmic scale to select an interval
    power10 = floor(log10(ideal_increment));
    base_power = 10^power10;
    
    % Standardized interval
    normalized = ideal_increment / base_power;
    
    % Extend the standard-interval array, Include finer decimal intervals
    % When the data range is very small, Use finer intervals
    standard_intervals = [0.001, 0.002, 0.005, ...
                          0.01, 0.02, 0.05, ...
                          0.1, 0.2, 0.5, ...
                          1, 2, 5, 10, 20, 50, 100, 200, 500, 1000];
    
    % Scale standard intervals to the current order of magnitude
    scaled_intervals = standard_intervals * 10^power10;
    
    % Find the standard interval nearest the target
    [~, idx] = min(abs(scaled_intervals - ideal_increment));
    increment = scaled_intervals(idx);
    
    % Ensure a nonzero interval
    if increment == 0
        increment = 0.1 * 10^power10;
    end
end

% Ensure the minimum increment is at least 0.001, To accommodate very small data ranges
if increment < 0.001
    increment = 0.001;
end

% Adjust y Axis maximum, Make it exactly a increment multiple
% Find the smallest value greater than or equal to y_max_plot that is a increment multiple
if increment > 0
    new_y_max = ceil(y_max_plot / increment) * increment;
else
    new_y_max = y_max_plot;
end

% If the new maximum is close to the current maximum, Adjust further
if (new_y_max - y_max_plot) < 0.3 * increment
    new_y_max = new_y_max + increment;
end

% For the minimum, Do not force adjustment to increment multiple
% Keep the minimum unchanged
new_y_min = y_min_plot;

% Adjust y axis range
set(ax, 'YLim', [new_y_min, new_y_max]);

% Generate evenly spaced ticks
% Redesign tick-generation logic, Ensure the entire y axis uses a uniform interval
ticks = [];

% From 0 generate ticks in both directions
if increment > 0
    % Ensure 0 Include ticks within the data range
    if ymin_data <= 0 && ymax_data >= 0
        % From 0 Generate ticks in the positive direction
        pos_ticks = 0:increment:new_y_max;
        % From 0 Generate ticks in the negative direction
        neg_ticks = 0:-increment:new_y_min;
    else
        % If 0 outside the data range, Start from the nearest 0 's increment multiple
        if new_y_min >= 0
            % Positive values only
            first_tick = ceil(new_y_min / increment) * increment;
            pos_ticks = first_tick:increment:new_y_max;
            neg_ticks = [];
        else
            % Negative values only
            last_tick = floor(new_y_max / increment) * increment;
            if new_y_max <= 0
                % Negative values only
                first_tick = floor(new_y_min / increment) * increment;
                neg_ticks = first_tick:increment:new_y_max;
                pos_ticks = [];
            else
                % Both positive and negative values, but 0 outside the range
                % From the last value less than or equal to new_y_min that is a increment multiple
                first_tick = floor(new_y_min / increment) * increment;
                % to the first value greater than or equal to new_y_max that is a increment multiple
                last_tick = ceil(new_y_max / increment) * increment;
                ticks = first_tick:increment:last_tick;
            end
        end
    end
    
    % Combine positive and negative ticks
    if isempty(ticks)
        if ~isempty(neg_ticks)
            ticks = [neg_ticks];
        end
        if ~isempty(pos_ticks)
            ticks = [ticks, pos_ticks];
        end
    end
    
    % Remove duplicate 0
    ticks = unique(round(ticks, 10));
    
    % Ensure ticks are within range
    ticks = ticks(ticks >= new_y_min & ticks <= new_y_max);
    
    % Ensure the maximum is included among ticks
    if ~ismember(new_y_max, ticks)
        ticks = [ticks, new_y_max];
    end
    
    % Check whether the minimum needs a label
    % Only when the minimum is a increment multiple, display its label
    if ~isempty(ticks)
        [min_tick, min_idx] = min(ticks);
        if abs(min_tick - new_y_min) < 1e-10
            % The minimum coincides with a tick, Check whether it is increment multiple
            remainder = mod(new_y_min, increment);
            if abs(remainder) > 1e-10 && abs(remainder - increment) > 1e-10
                % The minimum is not increment multiple, Mark the label as hidden
                show_min_label = false;
            else
                show_min_label = true;
            end
        else
            show_min_label = false;
        end
    else
        show_min_label = false;
    end
else
    ticks = linspace(new_y_min, new_y_max, desired_ticks);
    show_min_label = true;
end

% Sort
ticks = sort(ticks);

% Set ticks
set(ax, 'YTick', ticks);

% Compute the interval's decimal precision, For label formatting
if increment < 1
    decimal_places = max(0, -floor(log10(increment)));
else
    decimal_places = 0;
end

% Create tick labels
tick_labels = cell(1, length(ticks));
for i = 1:length(ticks)
    % Get current tick values
    tick_value = ticks(i);
    
    % Check whether to hide the minimum label
    if ~show_min_label && i == 1 && abs(tick_value - new_y_min) < 1e-10
        % The minimum does not follow uniform spacing, Hide label
        tick_labels{i} = '';
        continue;
    end
    
    % Check whether a tick is too close to its predecessor
    if i > 1 && abs(tick_value - ticks(i-1)) < 1e-6
        % If too close to the previous tick, Hide the label to avoid overlap
        tick_labels{i} = '';
        continue;
    end
    
    % Check whether a tick is an integer
    if abs(tick_value - round(tick_value)) < 1e-10
        % Integer, Use integer formatting
        tick_labels{i} = sprintf('%d', round(tick_value));
    else
        % Not an integer, Format with the specified decimal precision
        if decimal_places > 0
            % Format with the specified decimal precision
            format_str = sprintf('%%.%df', decimal_places);
            label_str = sprintf(format_str, tick_value);
            
            % Remove trailing 0
            label_str = regexprep(label_str, '0+$', '');
            % If a decimal point remains at the end, remove it as well
            label_str = regexprep(label_str, '\\.$', '');
        else
            % No decimal places, Display the integer directly
            label_str = sprintf('%d', round(tick_value));
        end
        
        tick_labels{i} = label_str;
    end
    
    % Fix negative-zero formatting
    if strcmp(tick_labels{i}, '-0')
        tick_labels{i} = '0';
    end
end

% Fix label formatting
for i = 1:length(tick_labels)
    % Check for multiple decimal points
    if ~isempty(tick_labels{i}) && contains(tick_labels{i}, '..')
        % If multiple decimal points occur, reformat
        if decimal_places > 0
            format_str = sprintf('%%.%df', decimal_places);
            tick_labels{i} = sprintf(format_str, ticks(i));
        else
            tick_labels{i} = sprintf('%d', round(ticks(i)));
        end
    end
    
    % Check whether the label equals its predecessor
    if i > 1 && ~isempty(tick_labels{i}) && strcmp(tick_labels{i}, tick_labels{i-1})
        % If equal to the previous label, Clear this label
        tick_labels{i} = '';
    end
end

% Set tick labels
set(ax, 'YTickLabel', tick_labels);
end

function ranovatbl = p2Ride2Ranova(xx1)
% Perform repeated-measures ANOVA
xx1_tmp = array2table(xx1, 'VariableNames', {'s1', 's2', 'p1', 'p2'});
rm = fitrm(xx1_tmp, 's1,s2,p1,p2 ~ 1', ...
           'WithinDesign', table([1 -1 1 -1]', [1 1 -1 -1]', ...
           'VariableNames', {'lambda'; 'drug'}));
ranovatbl = ranova(rm, 'WithinModel', 'lambda*drug');
end




function runDir = latest_run_local(parentDir)
    candidates = dir(parentDir);
    candidates = candidates([candidates.isdir] & ~ismember({candidates.name},{'.','..'}));
    assert(~isempty(candidates),'No fitted runs under %s. Run the corresponding fitting entry point first.',parentDir);
    [~,index] = max([candidates.datenum]);
    runDir = fullfile(candidates(index).folder,candidates(index).name);
end

function [os,ps,on,pn] = read_repeats_local(runDir)
    names = {'ots_data','pls_data','otn_data','pln_data'};
    allData = cell(1,4);
    for c = 1:4
        files = dir(fullfile(runDir,names{c},'repeat_*','results_all.mat'));
        assert(~isempty(files),'Missing repeated fits for %s under %s.',names{c},runDir);
        repeats = cell(numel(files),1);
        for r = 1:numel(files)
            fit = load(fullfile(files(r).folder,files(r).name),'results');
            repeats{r} = [cat(1,fit.results.optParams),cat(1,fit.results.bic),cat(1,fit.results.nll)];
            assert(size(repeats{r},1)==85,'Expected 85 completed participant fits.');
        end
        allData{c} = cell(85,1);
        for participant = 1:85
            allData{c}{participant} = cell2mat(cellfun(@(x)x(participant,:),repeats,'UniformOutput',false));
        end
    end
    os=allData{1};ps=allData{2};on=allData{3};pn=allData{4};
end

function n = count_repeats_local(runDir)
    files = dir(fullfile(runDir,'ots_data','*_*','repeat_*','results_all.mat'));
    assert(~isempty(files),'No completed comparison fits under %s.',runDir);
    n = max(cellfun(@(x)str2double(regexp(x,'(?<=repeat_)\d+$','match','once')), ...
        {files.folder}));
end
