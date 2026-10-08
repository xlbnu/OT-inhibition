function sig = plot_column_significance(x, data, varargin)
% PLOT_COLUMN_SIGNIFICANCE
%
% For data For each column, perform a one-sample t test, and mark significant time points on the current plot.
%
% Input: 
%   x    : 1 x n_time time axis
%   data : n_subject x n_time Data matrix
%
% Default behavior: 
%   - Each column ttest(data(:,i), 0)
%   - Connect adjacent significant points into a horizontal line
%   - Default significance-line position uses mean(data) ± SEM plot range
%
% Common examples: 
%   sig = plot_column_significance(vchan_time, beta_data, ...
%       'Color', [0.2 0.2 1], ...
%       'TimeWindow', [-0.8 0.5], ...
%       'Location', 'top');
%
% Optional parameters: 
%   'Alpha'           : Default 0.05
%   'Tail'            : 'both'/'right'/'left', Default 'both'
%   'Correction'      : 'none'/'fdr', Default 'none'
%   'Color'           : Significance-line color, Default [0.1 0.1 0.1]
%   'LineWidth'       : Significance-line width, Default 3
%   'TimeWindow'      : Mark significant points only within this window, All by default
%   'Location'        : 'top'/'bottom', Default 'top'
%   'YOffsetRatio'    : Distance between significance line and data boundary, Default 0.08
%   'YReference'      : 'summary'/'data'/'axis'/'value', Default 'summary'
%   'YValue'          : When YReference='value' manually specify y coordinates
%   'ConnectAdjacent' : true/false, Default true.Whether to connect adjacent significant points
%   'LineLength'      : Line length for an isolated significant point, Default 0.9 * median(diff(x))
%   'DeleteOld'       : true/false, Default false.Delete old lines with the same Tag significance marker
%   'Tag'             : Significance-line objects tag, Default 'column_significance_line'

p = inputParser;
p.addRequired('x', @(v) isnumeric(v) && isvector(v));
p.addRequired('data', @(v) isnumeric(v) && ismatrix(v));
p.addParameter('Alpha', 0.05, @(v) isnumeric(v) && isscalar(v));
p.addParameter('Tail', 'both', @(s) any(strcmpi(s, {'both','right','left'})));
p.addParameter('Correction', 'none', @(s) any(strcmpi(s, {'none','fdr'})));
p.addParameter('Color', [0.1 0.1 0.1], @(v) isnumeric(v) && numel(v) == 3);
p.addParameter('LineWidth', 3, @(v) isnumeric(v) && isscalar(v));
p.addParameter('TimeWindow', [], @(v) isempty(v) || (isnumeric(v) && numel(v) == 2));
p.addParameter('Location', 'top', @(s) any(strcmpi(s, {'top','bottom'})));
p.addParameter('YOffsetRatio', 0.08, @(v) isnumeric(v) && isscalar(v));
p.addParameter('YReference', 'summary', @(s) any(strcmpi(s, {'summary','data','axis','value'})));
p.addParameter('YValue', [], @(v) isempty(v) || (isnumeric(v) && isscalar(v)));
p.addParameter('ConnectAdjacent', true, @(v) islogical(v) && isscalar(v));
p.addParameter('LineLength', [], @(v) isempty(v) || (isnumeric(v) && isscalar(v)));
p.addParameter('AutoExpandYLim', true, @(v) islogical(v) && isscalar(v));
p.addParameter('DeleteOld', false, @(v) islogical(v) && isscalar(v));
p.addParameter('Tag', 'column_significance_line', @(s) ischar(s) || isstring(s));
p.parse(x, data, varargin{:});
opt = p.Results;

x = x(:).';

if size(data, 2) ~= numel(x) && size(data, 1) == numel(x)
    data = data.';
end
assert(size(data, 2) == numel(x), ...
    'data 必须是 subject x time，且列数等于 length(x)。');

n_time = numel(x);
pvals = nan(1, n_time);
tvals = nan(1, n_time);

if isempty(opt.TimeWindow)
    test_idx = true(1, n_time);
else
    test_idx = x >= opt.TimeWindow(1) & x <= opt.TimeWindow(2);
end

for i = find(test_idx)
    y = data(:, i);
    y = y(isfinite(y));
    if numel(y) >= 3 && std(y) > 0
        [~, p_this, ~, stat] = ttest(y, 0, ...
            'Alpha', opt.Alpha, 'Tail', opt.Tail);
        pvals(i) = p_this;
        tvals(i) = stat.tstat;
    end
end

qvals = nan(1, n_time);
h = false(1, n_time);

switch lower(opt.Correction)
    case 'none'
        h(test_idx) = pvals(test_idx) < opt.Alpha;
    case 'fdr'
        [h_tmp, ~, q_tmp] = simple_bh_fdr_local(pvals(test_idx), opt.Alpha);
        h(test_idx) = h_tmp;
        idx_tmp = find(test_idx);
        qvals(idx_tmp) = q_tmp;
end

ax = gca;
if opt.DeleteOld
    delete(findobj(ax, 'Tag', char(opt.Tag)));
end

yl_old = ylim(ax);
y_level = determine_y_level(data, test_idx, yl_old, opt);

if opt.AutoExpandYLim
    y_range = diff(yl_old);
    if y_range <= 0 || ~isfinite(y_range)
        y_range = 1;
    end
    switch lower(opt.Location)
        case 'top'
            if y_level >= yl_old(2)
                ylim(ax, [yl_old(1), y_level + 0.05*y_range]);
            end
        case 'bottom'
            if y_level <= yl_old(1)
                ylim(ax, [y_level - 0.05*y_range, yl_old(2)]);
            end
    end
end

dx = median(diff(x), 'omitnan');
if isempty(dx) || ~isfinite(dx) || dx <= 0
    dx = 1;
end

if isempty(opt.LineLength)
    line_len = 0.9 * dx;
else
    line_len = opt.LineLength;
end

hold_state = ishold(ax);
hold(ax, 'on');

sig_idx = find(h);
if ~isempty(sig_idx)
    if opt.ConnectAdjacent
        runs = contiguous_runs(sig_idx);
        for r = 1:size(runs, 1)
            idx1 = runs(r, 1);
            idx2 = runs(r, 2);

            x1 = x(idx1) - line_len/2;
            x2 = x(idx2) + line_len/2;

            plot(ax, [x1 x2], [y_level y_level], '-', ...
                'Color', opt.Color, ...
                'LineWidth', opt.LineWidth, ...
                'Tag', char(opt.Tag));
        end
    else
        for ii = 1:numel(sig_idx)
            i = sig_idx(ii);
            plot(ax, [x(i)-line_len/2, x(i)+line_len/2], ...
                [y_level, y_level], '-', ...
                'Color', opt.Color, ...
                'LineWidth', opt.LineWidth, ...
                'Tag', char(opt.Tag));
        end
    end
end

if ~hold_state
    hold(ax, 'off');
end

sig = struct();
sig.h = h;
sig.p = pvals;
sig.t = tvals;
sig.q = qvals;
sig.alpha = opt.Alpha;
sig.tail = opt.Tail;
sig.correction = opt.Correction;
sig.test_idx = test_idx;
sig.y_level = y_level;
sig.connect_adjacent = opt.ConnectAdjacent;

end


function y_level = determine_y_level(data, test_idx, yl_old, opt)

switch lower(opt.YReference)
    case 'summary'
        n_eff = sum(isfinite(data), 1);
        m = mean(data, 1, 'omitnan');
        se = std(data, 0, 1, 'omitnan') ./ sqrt(max(n_eff, 1));
        upper_y = m + se;
        lower_y = m - se;
        use_idx = test_idx & isfinite(upper_y) & isfinite(lower_y);
        if ~any(use_idx)
            use_idx = isfinite(upper_y) & isfinite(lower_y);
        end
        ref_max = max(upper_y(use_idx));
        ref_min = min(lower_y(use_idx));

    case 'data'
        d = data(:, test_idx);
        d = d(isfinite(d));
        if isempty(d)
            ref_min = yl_old(1);
            ref_max = yl_old(2);
        else
            ref_min = min(d);
            ref_max = max(d);
        end

    case 'axis'
        ref_min = yl_old(1);
        ref_max = yl_old(2);

    case 'value'
        assert(~isempty(opt.YValue), ...
            'YReference=''value'' 时必须提供 YValue。');
        y_level = opt.YValue;
        return;
end

ref_range = ref_max - ref_min;
if ~isfinite(ref_range) || ref_range <= 0
    ref_range = diff(yl_old);
end
if ~isfinite(ref_range) || ref_range <= 0
    ref_range = 1;
end

switch lower(opt.Location)
    case 'top'
        y_level = ref_max + opt.YOffsetRatio * ref_range;
    case 'bottom'
        y_level = ref_min - opt.YOffsetRatio * ref_range;
end

end


function runs = contiguous_runs(sig_idx)
% Return rows [start_index, end_index] in original time-index coordinates.

if isempty(sig_idx)
    runs = zeros(0, 2);
    return;
end

breaks = find(diff(sig_idx) > 1);
run_starts = [sig_idx(1), sig_idx(breaks + 1)];
run_ends = [sig_idx(breaks), sig_idx(end)];
runs = [run_starts(:), run_ends(:)];

end


function [h, crit_p, adj_p] = simple_bh_fdr_local(pvals, q)

p = pvals(:);
valid = isfinite(p);
pv = p(valid);

h = false(size(p));
adj_p = nan(size(p));
crit_p = 0;

if isempty(pv)
    h = reshape(h, size(pvals));
    adj_p = reshape(adj_p, size(pvals));
    return;
end

[ps, sort_idx] = sort(pv, 'ascend');
m = numel(ps);
thresh = ((1:m)' / m) * q;
below = ps <= thresh;

h_valid = false(size(pv));
if any(below)
    k = find(below, 1, 'last');
    crit_p = ps(k);
    h_valid = pv <= crit_p;
end

adj_sorted = ps .* m ./ (1:m)';
adj_sorted = flipud(cummin(flipud(adj_sorted)));
adj_sorted(adj_sorted > 1) = 1;

adj_valid = nan(size(pv));
adj_valid(sort_idx) = adj_sorted;

h(valid) = h_valid;
adj_p(valid) = adj_valid;

h = reshape(h, size(pvals));
adj_p = reshape(adj_p, size(pvals));

end

