function x_den = figureTmp1_nonGroup(xx, colors, user_cfg)
% figureTmp1_nonGroup  Draw non-group raincloud-style distributions.
%
%   x_den = figureTmp1_nonGroup(xx)
%   x_den = figureTmp1_nonGroup(xx, colors)
%   x_den = figureTmp1_nonGroup(xx, colors, user_cfg)
%
% Input
%   xx       : m-by-n numeric matrix. Each COLUMN is one distribution and
%              each row is one observation (e.g. 19-by-4 -> 4 plots).
%              NaN/Inf values are ignored separately within each column.
%   colors   : n-by-3 RGB matrix in [0,1]. If omitted, the first four
%              distributions use the default colors below. For n > 4,
%              the four colors are repeated cyclically.
%   user_cfg : optional structure. Its fields overwrite cfg below.
%
% Output
%   x_den    : x coordinates of the density baselines. They can be used
%              later to place labels or other annotations.
%
% Example
%   xx = randn(19,4) + (0:3)*0.4;
%   figure('Color','w');
%   x_den = figureTmp1_nonGroup(xx);
%   set(gca,'XTickLabel',{'A','B','C','D'});
%   ylabel('Value');

if nargin < 2 || isempty(colors)
    colors = [226 106  83; ... % primary red
               88 118 227; ... % primary blue
              222 190 184; ... % light red
              157 171 213] / 255; % light blue
end
if nargin < 3 || isempty(user_cfg)
    user_cfg = struct();
end

validateattributes(xx, {'numeric'}, {'2d','nonempty'}, mfilename, 'xx', 1);
if ~isreal(xx)
    error([mfilename ':ComplexData'], 'xx must contain real-valued data.');
end

n_plot = size(xx, 2);
if size(colors, 2) ~= 3 || isempty(colors) || any(~isfinite(colors(:))) || ...
        any(colors(:) < 0 | colors(:) > 1)
    error([mfilename ':InvalidColors'], ...
        'colors must be a nonempty k-by-3 RGB matrix with values in [0,1].');
end
if size(colors, 1) < n_plot
    colors = colors(mod(0:n_plot-1, size(colors,1)) + 1, :);
else
    colors = colors(1:n_plot, :);
end

% =========================================================
% 1. Layout and appearance configuration
% =========================================================
cfg.scale        = 1.0;  % global scale
cfg.start_x      = 2.0;  % first density baseline
cfg.w_density    = 2.5;  % maximum leftward density width
cfg.w_box        = 0.75; % box width
cfg.w_jitter     = 0.40; % scatter jitter width
cfg.pad_den_box  = 0.20; % density-to-box gap
cfg.pad_box_scat = 0.30; % box-to-scatter gap
cfg.pad_group    = 3.50; % gap between adjacent complete plots

cfg.num_density_points = 500;
cfg.marker             = 'o';
cfg.marker_size        = 105;
cfg.marker_alpha       = 0.70;
cfg.box_alpha          = 0.90;
cfg.line_width         = 2.0;
cfg.x_margin           = 0.50;
cfg.y_padding          = [0.05 0.08]; % [bottom top], fraction of data range
cfg.auto_y_ticks       = true;

cfg_fields = fieldnames(user_cfg);
for k = 1:numel(cfg_fields)
    cfg.(cfg_fields{k}) = user_cfg.(cfg_fields{k});
end

validateattributes(cfg.scale, {'numeric'}, {'scalar','real','positive','finite'});
validateattributes(cfg.y_padding, {'numeric'}, ...
    {'vector','numel',2,'real','nonnegative','finite'});

S = cfg.scale;
cfg.w_density    = cfg.w_density    * S;
cfg.w_box        = cfg.w_box        * S;
cfg.w_jitter     = cfg.w_jitter     * S;
cfg.pad_den_box  = cfg.pad_den_box  * S;
cfg.pad_box_scat = cfg.pad_box_scat * S;
cfg.pad_group    = cfg.pad_group    * S;
marker_size      = cfg.marker_size  * S^2;
line_width       = max(0.2, cfg.line_width * S);

% =========================================================
% 2. Equally spaced x coordinates: density -> box -> scatter
% =========================================================
x_den  = zeros(1, n_plot);
x_box  = zeros(1, n_plot);
x_scat = zeros(1, n_plot);
current_x = cfg.start_x;

for g = 1:n_plot
    x_den(g)  = current_x;
    x_box(g)  = x_den(g) + cfg.pad_den_box + cfg.w_box/2;
    x_scat(g) = x_box(g) + cfg.w_box/2 + ...
        cfg.pad_box_scat + cfg.w_jitter/2;
    current_x = x_scat(g) + cfg.w_jitter/2 + cfg.pad_group;
end

ax = gca;
hold_state = ishold(ax);
hold(ax, 'on');

finite_data = xx(isfinite(xx));
if isempty(finite_data)
    error([mfilename ':NoFiniteData'], ...
        'xx must contain at least one finite observation.');
end
data_min = min(finite_data);
data_max = max(finite_data);
data_range = data_max - data_min;
if data_range == 0
    data_range = max(1, abs(data_max));
end

% =========================================================
% 3. Left-facing density distributions
% =========================================================
for g = 1:n_plot
    values = xx(:, g);
    values = values(isfinite(values));
    if isempty(values)
        continue;
    end

    if numel(values) >= 2 && max(values) > min(values)
        [f, yi] = ksdensity(values, 'NumPoints', cfg.num_density_points);
        keep = yi >= min(values) & yi <= max(values);
        f = f(keep);
        yi = yi(keep);
        if ~isempty(f) && max(f) > 0
            f = f ./ max(f) .* cfg.w_density;
            fill(ax, x_den(g) - [0 f 0], ...
                [min(values) yi max(values)], colors(g,:), ...
                'EdgeColor', colors(g,:), 'FaceAlpha', 1, ...
                'LineWidth', line_width);
        end
    else
        % A narrow symmetric glyph keeps a constant/single-point row visible.
        y0 = values(1);
        dy = 0.015 * data_range;
        fill(ax, x_den(g) - [0 cfg.w_density 0], ...
            [y0-dy y0 y0+dy], colors(g,:), ...
            'EdgeColor', colors(g,:), 'FaceAlpha', 1, ...
            'LineWidth', line_width);
    end
end

% =========================================================
% 4. Boxplots (one column per input row)
% =========================================================
% Each input column is already one distribution for boxplot.
% boxplot(ax, xx, 'Positions', x_box, 'Colors', [1 1 1], ...
%     'Widths', cfg.w_box, 'Symbol', '');
% box_objects = findobj(ax, 'Tag', 'Box');
% for k = 1:numel(box_objects)
%     xd = get(box_objects(k), 'XData');
%     [~, g] = min(abs(mean(xd) - x_box));
%     patch(ax, xd, get(box_objects(k), 'YData'), colors(g,:), ...
%         'EdgeColor', 'none', 'FaceAlpha', cfg.box_alpha);
% end
% box_lines = boxplot(ax, xx, 'Positions', x_box, 'Colors', [0 0 0], ...
%     'Widths', cfg.w_box, 'Symbol', '');
% set(box_lines, 'LineWidth', line_width);

% =========================================================
% 5. Adaptively sized box plots
% =========================================================
boxplot(xx, 'Positions', x_box, 'Color', [1 1 1], 'Widths', cfg.w_box, 'Symbol', '');
boxObj = findobj(gca, 'Tag', 'Box'); 
for i = 1:length(boxObj)
    patch(get(boxObj(i), 'XData'), get(boxObj(i), 'YData'), colors(5-i,:), 'FaceAlpha', 0.9);
end
bx = boxplot(xx, 'Positions', x_box, 'color', [0 0 0], 'Widths', cfg.w_box ,'Symbol', '');
set(bx,'LineWidth', line_width);
try RemoveBoxplotlines(gca); catch; end 

% =========================================================
% 5. Independent scatter points (no pairwise connecting lines)
% =========================================================
for g = 1:n_plot
    values = xx(:, g);
    valid = isfinite(values);
    n_valid = nnz(valid);
    if n_valid == 0
        continue;
    end
    jitter_x = x_scat(g) + (rand(n_valid, 1) - 0.5) .* cfg.w_jitter;
    scatter(ax, jitter_x, values(valid), marker_size, cfg.marker, ...
        'MarkerFaceColor', colors(g,:), 'MarkerEdgeColor', [1 1 1], ...
        'MarkerFaceAlpha', cfg.marker_alpha);
end

% =========================================================
% 6. Axes: compact automatic y limits and clean rendering
% =========================================================
y_min_plot = data_min - cfg.y_padding(1) * data_range;
y_max_plot = data_max + cfg.y_padding(2) * data_range;
if cfg.auto_y_ticks
    set_nice_y_ticks(ax, y_min_plot, y_max_plot);
else
    ylim(ax, [y_min_plot y_max_plot]);
end

left_edge  = x_den(1) - cfg.w_density - cfg.x_margin*S;
right_edge = x_scat(end) + cfg.w_jitter/2 + cfg.x_margin*S;
xlim(ax, [left_edge right_edge]);
box(ax, 'off');
grid(ax, 'off');
set(ax, 'FontName', 'Arial', 'TickDir', 'in', ...
    'TickLabelInterpreter', 'none', 'LineWidth', 3);
set(ancestor(ax, 'figure'), 'Renderer', 'painters');

if ~hold_state
    hold(ax, 'off');
end
end

function set_nice_y_ticks(ax, y_min, y_max)
% Select 4-7 readable ticks while retaining the requested data padding.
span = y_max - y_min;
if ~(isfinite(span) && span > 0)
    y_min = -1;
    y_max = 1;
    span = 2;
end

raw_step = span / 5;
power10 = 10^floor(log10(raw_step));
candidates = [1 2 2.5 5 10] * power10;
best_score = inf;
best_ticks = [];

for step = candidates
    lo = floor(y_min / step) * step;
    hi = ceil(y_max / step) * step;
    ticks = lo:step:(hi + step*0.01);
    n_tick = numel(ticks);
    if n_tick < 4 || n_tick > 7
        continue;
    end
    score = abs(n_tick - 6) + 0.2 * ((hi-lo)-span) / span;
    if score < best_score
        best_score = score;
        best_ticks = ticks;
    end
end

if isempty(best_ticks)
    step = candidates(2);
    best_ticks = floor(y_min/step)*step:step:ceil(y_max/step)*step;
end
best_ticks(abs(best_ticks) < eps(max(1,max(abs(best_ticks))))*10) = 0;
ylim(ax, [best_ticks(1) best_ticks(end)]);
yticks(ax, best_ticks);
end
