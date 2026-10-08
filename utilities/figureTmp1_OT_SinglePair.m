function x_den = figureTmp1_OT_SinglePair(xx1, x_num, ot_color, user_cfg)
% figureTmp1_OT_SinglePair - Plot a single paired group distribution (N x 2)
% including: Density plot (Outward-facing on both sides) + Boxplot + Connected jittered scatter plot

if nargin < 2, x_num = 1; end
if nargin < 3 || isempty(ot_color)
    % Use only the first two original colors by default
    ot_color = [226 106 83; 242 199 178]./255; 
end

% Use the first two colors to avoid out-of-bounds indexing
ot_color = ot_color(1:2, :);

% Get data length, For adaptive layout and shape checks
data_len = size(xx1, 1);

% =========================================================
% 1. Flexible layout configuration panel
% =========================================================
cfg.scale         = 1.0;   % Global scale factor
cfg.start_x       = 2.0;   % Base starting position X coordinates
cfg.w_density     = 2.5;   % Fill Outward extension width
cfg.w_box         = 0.75;  % Boxplot width
cfg.w_jitter      = 0.4;   % Scatter Jitter-width range
cfg.pad_den_box   = 0.2;   % Density Up to Boxplot distances
cfg.pad_box_scat  = 0.3;   % Boxplot Up to Scatter spacing
cfg.pad_scat_scat = 0.2;   % 3 wave Scatter spacing between (85 Dedicated)
cfg.pad_pair      = 2.0;   % Gap between connecting lines within pairs

if nargin == 4
    fields = fieldnames(user_cfg);
    for i = 1:numel(fields), cfg.(fields{i}) = user_cfg.(fields{i}); end
end

% Apply scale multipliers
S = cfg.scale;
cfg.w_density     = cfg.w_density * S;
cfg.w_box         = cfg.w_box * S;
cfg.w_jitter      = cfg.w_jitter * S;
cfg.pad_den_box   = cfg.pad_den_box * S;
cfg.pad_box_scat  = cfg.pad_box_scat * S;
cfg.pad_scat_scat = cfg.pad_scat_scat * S;
cfg.pad_pair      = cfg.pad_pair * S;

sz_cir = 55 * (S^2); 
sz_squ = 60 * (S^2);
sz_tri = 58 * (S^2);
lw_line = max(0.1, 0.4 * S); 
lw_box  = max(0.2, 2.0 * S);
sig_textSize = 33;
lw_line_text = 3;

hold on;

% =========================================================
% 2. Dynamically compute all elements' X coordinate centers (Single-group mirror symmetry)
% =========================================================
x_den = zeros(1, 2);
x_box = zeros(1, 2);
x_scat_cir = zeros(1, 2);
x_scat_squ = zeros(1, 2);
x_scat_tri = zeros(1, 2);

c1 = 1; c2 = 2;

% Compute the left side (variables 1) coordinates
x_den(c1) = cfg.start_x; 
x_box(c1) = x_den(c1) + cfg.pad_den_box + cfg.w_box/2; 

if data_len == 85
    x_scat_cir(c1) = x_box(c1) + cfg.w_box/2 + cfg.pad_box_scat + cfg.w_jitter/2;
    x_scat_squ(c1) = x_scat_cir(c1) + cfg.w_jitter/2 + cfg.pad_scat_scat + cfg.w_jitter/2;
    x_scat_tri(c1) = x_scat_squ(c1) + cfg.w_jitter/2 + cfg.pad_scat_scat + cfg.w_jitter/2;
    right_edge_c1  = x_scat_tri(c1) + cfg.w_jitter/2;
else
    x_single = x_box(c1) + cfg.w_box/2 + cfg.pad_box_scat + cfg.w_jitter/2;
    x_scat_cir(c1) = x_single;
    x_scat_squ(c1) = x_single;
    x_scat_tri(c1) = x_single;
    right_edge_c1  = x_single + cfg.w_jitter/2;
end

left_edge_c2 = right_edge_c1 + cfg.pad_pair;

% Compute the right side (variables 2) coordinates, Mirror-symmetric arrangement
if data_len == 85
    x_scat_tri(c2) = left_edge_c2 + cfg.w_jitter/2; 
    x_scat_squ(c2) = x_scat_tri(c2) + cfg.w_jitter/2 + cfg.pad_scat_scat + cfg.w_jitter/2;
    x_scat_cir(c2) = x_scat_squ(c2) + cfg.w_jitter/2 + cfg.pad_scat_scat + cfg.w_jitter/2;
    x_box(c2) = x_scat_cir(c2) + cfg.w_jitter/2 + cfg.pad_box_scat + cfg.w_box/2;
else
    x_single = left_edge_c2 + cfg.w_jitter/2;
    x_scat_tri(c2) = x_single;
    x_scat_squ(c2) = x_single;
    x_scat_cir(c2) = x_single;
    x_box(c2) = x_single + cfg.w_jitter/2 + cfg.pad_box_scat + cfg.w_box/2;
end

x_den(c2) = x_box(c2) + cfg.w_box/2 + cfg.pad_den_box;

% =========================================================
% 3. Plot kernel-density estimates (Fill) - Left side extends leftward, Right side extends rightward
% =========================================================
% Left side (c1)
x_tmp = xx1(:, c1);
[f, xi] = ksdensity(x_tmp, 'NumPoints', 500);
f = f ./ max(f) .* cfg.w_density; 
ql0 = max(x_tmp) .* 1.01; qu0 = min(x_tmp) .* 0.98;
xu = find(xi > ql0, 1); xl = find(xi < qu0, 1, 'last');
f(xi > ql0 | xi < qu0) = 0;
xi(xu) = xi(xu-1); xi(xl) = xi(xl+1);
xi = xi(xl:xu); f = f(xl:xu);
fill(-f + x_den(c1), xi, ot_color(c1,:), 'EdgeColor', [0 0 0], 'FaceAlpha', 1,'LineWidth',lw_box);

% Right side (c2)
x_tmp = xx1(:, c2);
[f, xi] = ksdensity(x_tmp, 'NumPoints', 500);
f = f ./ max(f) .* cfg.w_density;
ql0 = max(x_tmp) .* 1.01; qu0 = min(x_tmp) .* 0.98;
xu = find(xi > ql0, 1); xl = find(xi < qu0, 1, 'last');
f(xi > ql0 | xi < qu0) = 0;
xi(xu) = xi(xu-1); xi(xl) = xi(xl+1);
xi = xi(xl:xu); f = f(xl:xu);
fill(f + x_den(c2), xi, ot_color(c2,:), 'EdgeColor', [0 0 0], 'FaceAlpha', 1,'LineWidth',lw_box);

% =========================================================
% 4. Generate scatter points and connecting lines
% =========================================================
x_pos1 = zeros(data_len, 2);
for c = 1:2
    if data_len == 85
        x_pos1(1:37, c)  = (rand(37, 1) - 0.5) * cfg.w_jitter + x_scat_cir(c);
        x_pos1(38:56, c) = (rand(19, 1) - 0.5) * cfg.w_jitter + x_scat_squ(c);
        x_pos1(57:85, c) = (rand(29, 1) - 0.5) * cfg.w_jitter + x_scat_tri(c);
    elseif data_len == 37
        x_pos1(1:37, c)  = (rand(37, 1) - 0.5) * cfg.w_jitter + x_scat_cir(c);
    elseif data_len == 19
        x_pos1(1:19, c)  = (rand(19, 1) - 0.5) * cfg.w_jitter + x_scat_squ(c);
    elseif data_len == 29
        x_pos1(1:29, c)  = (rand(29, 1) - 0.5) * cfg.w_jitter + x_scat_tri(c);
    else
        % Handle general cases
        x_pos1(:, c) = (rand(data_len, 1) - 0.5) * cfg.w_jitter + x_scat_cir(c);
    end
end

% Plot paired connecting lines
pl1 = plot(x_pos1', xx1', 'Color', [0.7 0.7 0.7], 'LineWidth', lw_line);
for i = 1:length(pl1)
    pl1(i).Color(4) = 0.5;
end

% Plot foreground scatter points
if data_len == 85
    sc1 = scatter(x_pos1(1:37,:), xx1(1:37, :), sz_cir, 'o', 'MarkerFaceColor', 'flat', 'MarkerEdgeColor', [1 1 1], 'MarkerFaceAlpha', 0.7);
    sc2 = scatter(x_pos1(38:56,:), xx1(38:56, :), sz_squ, 'square', 'MarkerFaceColor', 'flat', 'MarkerEdgeColor', [1 1 1], 'MarkerFaceAlpha', 0.7);
    sc3 = scatter(x_pos1(57:85,:), xx1(57:85, :), sz_tri, '^', 'MarkerFaceColor', 'flat', 'MarkerEdgeColor', [1 1 1], 'MarkerFaceAlpha', 0.7);
    for i = 1:2
        sc1(i).MarkerFaceColor = ot_color(i,:);
        sc2(i).MarkerFaceColor = ot_color(i,:);
        sc3(i).MarkerFaceColor = ot_color(i,:);
    end
elseif data_len == 37
    sc1 = scatter(x_pos1, xx1, sz_cir+50, 'o', 'MarkerFaceColor', 'flat', 'MarkerEdgeColor', [1 1 1], 'MarkerFaceAlpha', 0.7);
    for i = 1:2, sc1(i).MarkerFaceColor = ot_color(i,:); end
elseif data_len == 19
    sc2 = scatter(x_pos1, xx1, sz_squ+50, 'square', 'MarkerFaceColor', 'flat', 'MarkerEdgeColor', [1 1 1], 'MarkerFaceAlpha', 0.7);
    for i = 1:2, sc2(i).MarkerFaceColor = ot_color(i,:); end
elseif data_len == 29
    sc3 = scatter(x_pos1, xx1, sz_tri+50, '^', 'MarkerFaceColor', 'flat', 'MarkerEdgeColor', [1 1 1], 'MarkerFaceAlpha', 0.7);
    for i = 1:2, sc3(i).MarkerFaceColor = ot_color(i,:); end
else
    % General plotting
    sc_other = scatter(x_pos1, xx1, sz_cir+50, 'o', 'MarkerFaceColor', 'flat', 'MarkerEdgeColor', [1 1 1], 'MarkerFaceAlpha', 0.7);
    for i = 1:2, sc_other(i).MarkerFaceColor = ot_color(i,:); end
end

% =========================================================
% 5. Dynamically adjusted boxplots
% =========================================================
boxplot(xx1, 'Positions', x_box, 'Color', [1 1 1], 'Widths', cfg.w_box, 'Symbol', '');
boxObj = findobj(gca, 'Tag', 'Box'); 

% Note: Handle Graphics The stack is typically reversed.If present 2 items Box, Index is 1 Represents the right side (c2), Index is 2 Represents the left side (c1)
for i = 1:length(boxObj)
    patch(get(boxObj(i), 'XData'), get(boxObj(i), 'YData'), ot_color(3-i,:), 'FaceAlpha', 0.9);
end
bx = boxplot(xx1, 'Positions', x_box, 'color', [0 0 0], 'Widths', cfg.w_box ,'Symbol', '');
set(bx,'LineWidth', lw_box);
try RemoveBoxplotlines(gca); catch; end 

% =========================================================
% 6. Significance calculation and Y adaptive axis limits
% =========================================================
[~, p] = ttest(xx1(:, 1), xx1(:, 2));
p_adj = p .* x_num;

if p_adj < 0.001, sig_text = '***';
elseif p_adj < 0.01, sig_text = '**';
elseif p_adj < 0.05, sig_text = '*';
elseif p_adj < 0.08, sig_text = ['p=',mat2str(round(p_adj,4))];
else, sig_text = 'n.s.';
end

% Compute the actual range of original data
ymax_data = max(xx1, [], 'all');
ymin_data = min(xx1, [], 'all');
y_range_data = ymax_data - ymin_data;
if y_range_data == 0, y_range_data = 1; end

% Provide 20%base space to find suitable ticks
y_min_plot = ymin_data - 0.05 * y_range_data;
y_max_plot = ymax_data + 0.20 * y_range_data;  

% Call the latest adaptive natural-number tick selector
[~, new_y_max, ~] = setCustomYTicks(gca, y_min_plot, y_max_plot, ymax_data);

% =========================================================
% ★ Flexible vertical layout (Flex-Layout)
% =========================================================
top_space = new_y_max - ymax_data; 
sig_level1_y    = ymax_data + 0.35 * top_space; 
sig_text_offset = 0.01 * top_space;             

x_pair = [x_den(1), x_den(2)];
plot(x_pair, [sig_level1_y, sig_level1_y], 'k-', 'LineWidth', lw_line_text);

if contains(sig_text, '**')
    sig_textSize_tmp = sig_textSize + 10;
    y_text_pos = sig_level1_y + sig_text_offset - 0.20 * top_space; 
elseif contains(sig_text, '*')
    sig_textSize_tmp = sig_textSize + 10;
    y_text_pos = sig_level1_y + sig_text_offset - 0.15 * top_space;         
else
    sig_textSize_tmp = sig_textSize;
    y_text_pos = sig_level1_y + sig_text_offset;
end

text(mean(x_pair), y_text_pos, sig_text, ...
     'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
     'FontSize', sig_textSize_tmp, 'FontName', 'Arial', 'Interpreter', 'none');

% =========================================================
% ★ 7. Protect text and low-level graphical rendering
% =========================================================
box off;
try gcaf1(gca); catch; end
grid off;

% Dynamically scale X axis to the width of a single-group plot
xlim([x_den(1) - cfg.w_density - 1 * S, x_den(2) + cfg.w_density + 1 * S]);

set(gcf, 'Renderer', 'painters');
set(gca, 'FontName', 'Arial');
set(gca, 'TickDir', 'out', 'TickLabelInterpreter', 'none'); 
set(gca, 'LineWidth', 3, 'TickDir', 'in'); 
end

%% ====== Y Axis tick-control component (Keep unchanged) ======
function [new_y_min, new_y_max, increment] = setCustomYTicks(ax, y_min_plot, y_max_plot, ymax_data)
    y_range = y_max_plot - y_min_plot;
    if y_range <= 0 || isnan(y_range) || isinf(y_range)
        y_min_plot = -1; y_max_plot = 1; ymax_data = 0.5;
        y_range = 2;
    end
    
    mag = 10^floor(log10(y_range / 5)); 
    candidate_ds = [];
    for m = [mag/10, mag, mag*10]
        for n = [1, 1.5, 2, 2.5, 3, 4, 5, 6, 7, 8, 9]
            candidate_ds(end+1) = n * m;
        end
    end
    
    best_score = inf;
    best_d = candidate_ds(1);
    best_max = y_max_plot;
    best_ticks = [];
    
    for i = 1:length(candidate_ds)
        d = candidate_ds(i);
        test_max = ceil(y_max_plot / d) * d;
        first_tick = ceil(y_min_plot / d) * d;
        
        if first_tick > test_max, continue; end
        
        ticks = first_tick : d : (test_max + d * 0.01);
        ticks = round(ticks ./ d) .* d;
        n_ticks = length(ticks);
        
        if n_ticks < 3 || n_ticks > 8, continue; end
        
        valid_ticks_flag = true;
        for t = ticks
            v = abs(t);
            if v < 1e-10 * d, continue; end 
            s = num2str(v, '%g');
            e_idx = find(s == 'e' | s == 'E');
            if ~isempty(e_idx), s = s(1:e_idx-1); end
            s = strrep(s, '.', '');
            s = regexprep(s, '^0+', ''); 
            s = regexprep(s, '0+$', ''); 
            
            if length(s) > 1
                if ~ismember(s(end), ['0','2','4','6','8','5'])
                    valid_ticks_flag = false;
                    break;
                end
            end
        end
        if ~valid_ticks_flag, continue; end
        
        stretch = (test_max - y_max_plot) / y_range; 
        tick_penalty = max(0, n_ticks - 6) * 0.05;   
        
        top_ratio = (test_max - ymax_data) / (test_max - y_min_plot);
        top_penalty = 0;
        if top_ratio > 0.25
            top_penalty = (top_ratio - 0.25) * 50; 
        end
        
        base_n = round(d / (10^floor(log10(d))));
        if ismember(base_n, [1, 2, 5]), nice_penalty = 0;
        elseif ismember(base_n, [4, 8]), nice_penalty = 0.01;
        else, nice_penalty = 0.02; end
        
        score = stretch + tick_penalty + top_penalty + nice_penalty;
        
        if score < best_score
            best_score = score;
            best_d = d;
            best_max = test_max;
            best_ticks = ticks;
        end
    end
    
    if isempty(best_ticks)
        best_d = mag * 2; 
        best_max = ceil(y_max_plot / best_d) * best_d;
        first_tick = ceil(y_min_plot / best_d) * best_d;
        best_ticks = first_tick : best_d : (best_max + best_d * 0.01);
        best_ticks = round(best_ticks ./ best_d) .* best_d;
    end
    
    increment = best_d;
    new_y_min = y_min_plot; 
    new_y_max = best_max;
    
    set(ax, 'YLim', [new_y_min, new_y_max]);
    set(ax, 'YTick', best_ticks);
    
    tick_labels = cell(1, length(best_ticks));
    for i = 1:length(best_ticks)
        if abs(best_ticks(i)) < 1e-10 * increment, best_ticks(i) = 0; end
        tick_labels{i} = num2str(best_ticks(i), '%g');
    end
    
    set(ax, 'YTickLabel', tick_labels);
end