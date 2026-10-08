function x_den=figureTmp1_OT(xx1, x_num, ot_color, user_cfg)
% figureTmp1 - Plot paired data distributions (Advanced layout: Strict control of margins + Adaptive significant-digit filtering)
if nargin < 2, x_num = 1; end
if nargin < 3 || isempty(ot_color)
    ot_color = [226 106 83;242 199 178; ...
                88 118 227;171 198 251]./255;
end

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
cfg.pad_scat_scat = 0.2;   % 3 wave Scatter spacing between
cfg.pad_pair      = 2.0;   % Paired gap
cfg.pad_group     = 5.0;   % Gap between major groups

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
cfg.pad_group     = cfg.pad_group * S;

sz_cir = 55 * (S^2); 
sz_squ = 60 * (S^2);
sz_tri = 58 * (S^2);
lw_line = max(0.1, 0.4 * S); 
lw_box  = max(0.2, 2.0 * S);
sig_textSize = 33;
lw_line_text = 3;

hold on;

% =========================================================
% 2. Dynamically compute all elements' X coordinate centers
% =========================================================
x_den = zeros(1, 4);
x_box = zeros(1, 4);
x_scat_cir = zeros(1, 4);
x_scat_squ = zeros(1, 4);
x_scat_tri = zeros(1, 4);
current_x = cfg.start_x;

for i = 1:2
    c1 = 2*i - 1; c2 = 2*i;
    
    x_den(c1) = current_x; 
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
    
    if i < 2
        current_x = x_den(c2) + cfg.w_density + cfg.pad_group;
    end
end

% =========================================================
% 3. Plot kernel-density estimates (Fill)
% =========================================================
for i = 1:2
    c1 = 2*i - 1; c2 = 2*i;
    
    x_tmp = xx1(:, c1);
    [f, xi] = ksdensity(x_tmp, 'NumPoints', 500);
    f = f ./ max(f) .* cfg.w_density; 
    ql0 = max(x_tmp) .* 1.01; qu0 = min(x_tmp) .* 0.98;
    xu = find(xi > ql0, 1); xl = find(xi < qu0, 1, 'last');
    f(xi > ql0 | xi < qu0) = 0;
    xi(xu) = xi(xu-1); xi(xl) = xi(xl+1);
    xi = xi(xl:xu); f = f(xl:xu);
    fill(-f + x_den(c1), xi, ot_color(c1,:), 'EdgeColor', ot_color(c1,:), 'FaceAlpha', 1,'LineWidth',lw_box);
    
    x_tmp = xx1(:, c2);
    [f, xi] = ksdensity(x_tmp, 'NumPoints', 500);
    f = f ./ max(f) .* cfg.w_density;
    ql0 = max(x_tmp) .* 1.01; qu0 = min(x_tmp) .* 0.98;
    xu = find(xi > ql0, 1); xl = find(xi < qu0, 1, 'last');
    f(xi > ql0 | xi < qu0) = 0;
    xi(xu) = xi(xu-1); xi(xl) = xi(xl+1);
    xi = xi(xl:xu); f = f(xl:xu);
    fill(f + x_den(c2), xi, ot_color(c2,:), 'EdgeColor', ot_color(c2,:), 'FaceAlpha', 1,'LineWidth',lw_box);
end

% =========================================================
% 4. Generate and plot scatter matrices
% =========================================================
x_pos1 = zeros(data_len, 4);
for c = 1:4
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
    end
end
pl1 = plot(x_pos1(:, 1:2)', xx1(:, [1 2])', 'Color', [0.7 0.7 0.7], 'LineWidth', lw_line);
pl2 = plot(x_pos1(:, 3:4)', xx1(:, [3 4])', 'Color', [0.7 0.7 0.7], 'LineWidth', lw_line);
for i = 1:length(pl1)
    pl1(i).Color(4) = 0.5; pl2(i).Color(4) = 0.5;
end
if data_len == 85
    sc1 = scatter(x_pos1(1:37,:), xx1(1:37, :), sz_cir, 'o', 'MarkerFaceColor', 'flat', 'MarkerEdgeColor', [1 1 1], 'MarkerFaceAlpha', 0.7);
    sc2 = scatter(x_pos1(38:56,:), xx1(38:56, :), sz_squ, 'square', 'MarkerFaceColor', 'flat', 'MarkerEdgeColor', [1 1 1], 'MarkerFaceAlpha', 0.7);
    sc3 = scatter(x_pos1(57:85,:), xx1(57:85, :), sz_tri, '^', 'MarkerFaceColor', 'flat', 'MarkerEdgeColor', [1 1 1], 'MarkerFaceAlpha', 0.7);
    for i = 1:4
        sc1(i).MarkerFaceColor = ot_color(i,:);
        sc2(i).MarkerFaceColor = ot_color(i,:);
        sc3(i).MarkerFaceColor = ot_color(i,:);
    end
elseif data_len == 37
    sc1 = scatter(x_pos1, xx1, sz_cir+50, 'o', 'MarkerFaceColor', 'flat', 'MarkerEdgeColor', [1 1 1], 'MarkerFaceAlpha', 0.7);
    for i = 1:4, sc1(i).MarkerFaceColor = ot_color(i,:); end
elseif data_len == 19
    sc2 = scatter(x_pos1, xx1, sz_squ+50, 'square', 'MarkerFaceColor', 'flat', 'MarkerEdgeColor', [1 1 1], 'MarkerFaceAlpha', 0.7);
    for i = 1:4, sc2(i).MarkerFaceColor = ot_color(i,:); end
elseif data_len == 29
    sc3 = scatter(x_pos1, xx1, sz_tri+50, '^', 'MarkerFaceColor', 'flat', 'MarkerEdgeColor', [1 1 1], 'MarkerFaceAlpha', 0.7);
    for i = 1:4, sc3(i).MarkerFaceColor = ot_color(i,:); end
end

% =========================================================
% 5. Dynamically adjusted boxplots
% =========================================================
boxplot(xx1, 'Positions', x_box, 'Color', [1 1 1], 'Widths', cfg.w_box, 'Symbol', '');
boxObj = findobj(gca, 'Tag', 'Box'); 
for i = 1:length(boxObj)
    patch(get(boxObj(i), 'XData'), get(boxObj(i), 'YData'), ot_color(5-i,:), 'FaceAlpha', 0.9);
end
bx = boxplot(xx1, 'Positions', x_box, 'color', [0 0 0], 'Widths', cfg.w_box ,'Symbol', '');
set(bx,'LineWidth', lw_box);
try RemoveBoxplotlines(gca); catch; end 

% =========================================================
% 6. Significance calculation and Y adaptive axis limits
% =========================================================
sig_texts = cell(1, 2);
for i = 1:2
    pair_idx = 2*i-1:2*i;
    [~, p] = ttest(xx1(:, pair_idx(1)), xx1(:, pair_idx(2)));
    p_adj = p.*x_num;
    if p_adj < 0.001, sig_texts{i} = '***';
    elseif p_adj < 0.01, sig_texts{i} = '**';
    elseif p_adj < 0.05, sig_texts{i} = '*';
    % elseif p_adj < 0.08, sig_texts{i} = ['p=',mat2str(round(p_adj,4))];
    else, sig_texts{i} = 'n.s.';
    end
    % sig_texts{i} = ['p=',mat2str(round(p_adj,4))];
end
ranovatbl = p2Ride2Ranova(xx1);
re_padj = ranovatbl.pValue(7).*x_num;
if re_padj < 0.001, between_sig = '***';
elseif re_padj < 0.01, between_sig = '**';
elseif re_padj < 0.05, between_sig = '*';
elseif re_padj < 0.08, between_sig = ['p=',mat2str(round(re_padj,4))];  
else, between_sig = 'n.s.';
end
% between_sig = ['p=',mat2str(round(re_padj,4))];
% Compute the actual range of original data
ymax_data = max(xx1, [], 'all');
ymin_data = min(xx1, [], 'all');
y_range_data = ymax_data - ymin_data;
if y_range_data == 0, y_range_data = 1; end

% Provide 20%base space to find suitable ticks
y_min_plot = ymin_data - 0.05 * y_range_data;
y_max_plot = ymax_data + 0.20 * y_range_data;  

% Call the latest adaptive natural-number tick selector (Includes significant-digit checks and reference-line control)
[~, new_y_max, ~] = setCustomYTicks(gca, y_min_plot, y_max_plot, ymax_data);

% =========================================================
% ★ Flexible vertical layout (Flex-Layout): Allocate adaptively to available space
% =========================================================
top_space = new_y_max - ymax_data; 
sig_level1_y    = ymax_data + 0.35 * top_space; 
sig_level2_y    = ymax_data + 0.8 * top_space; 
sig_text_offset = 0.01 * top_space;             

for i = 1:2
    c1 = 2*i-1; c2 = 2*i;
    x_pair = [x_den(c1), x_den(c2)];
    plot(x_pair, [sig_level1_y, sig_level1_y], 'k-', 'LineWidth', lw_line_text);
    
    if contains(sig_texts{i}, '**')
        sig_textSize_tmp = sig_textSize + 10;
        y_text_pos = sig_level1_y + sig_text_offset - 0.175 * top_space; 
    elseif contains(sig_texts{i}, '*')
        sig_textSize_tmp = sig_textSize + 10;
        y_text_pos = sig_level1_y + sig_text_offset - 0.175 * top_space;         
    else
        sig_textSize_tmp = sig_textSize;
        y_text_pos = sig_level1_y + sig_text_offset;
    end
    
    text(mean(x_pair), y_text_pos, sig_texts{i}, ...
         'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
         'FontSize', sig_textSize_tmp, 'FontName', 'Arial', 'Interpreter', 'none');
end

x_group1 = mean([x_den(1), x_den(2)]);
x_group2 = mean([x_den(3), x_den(4)]);
plot([x_group1, x_group2], [sig_level2_y, sig_level2_y], 'k-', 'LineWidth', lw_line_text);

if contains(between_sig, '**')
    sig_textSize_tmp = sig_textSize + 10;
    y_text_pos = sig_level2_y + sig_text_offset - 0.175 * top_space;
elseif contains(between_sig, '*')
    sig_textSize_tmp = sig_textSize + 10;
    y_text_pos = sig_level2_y + sig_text_offset - 0.175 * top_space;
else
    sig_textSize_tmp = sig_textSize;
    y_text_pos = sig_level2_y + sig_text_offset;
end
text(mean([x_group1, x_group2]), y_text_pos, between_sig, ...
     'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', ...
     'FontSize', sig_textSize_tmp, 'FontName', 'Arial', 'Interpreter', 'none');

% =========================================================
% ★ 7. Protect text and low-level graphical rendering
% =========================================================
box off;
try gcaf1(gca); catch; end
grid off;
xlim([cfg.start_x - cfg.pad_group.*2/3 - 0.5 * S, x_den(4) + cfg.pad_group.*2/3 + 0.5 * S]);
set(gcf, 'Renderer', 'painters');
set(gca, 'FontName', 'Arial');
set(gca, 'TickDir', 'out', 'TickLabelInterpreter', 'none'); 
set(gca, 'LineWidth', 3, 'TickDir', 'in'); 
end

function setAxesForPPT(ax, fontSize_pt, axesLineWidth_pt)
% Apply final axis styling and use U+2212 for negative tick labels.
% Positive tick labels are left unchanged.
    set(ax, 'FontName','Arial', 'FontSize',fontSize_pt, ...
        'LineWidth',axesLineWidth_pt, 'TickLabelInterpreter','none');
    set(ax.XLabel, 'FontName','Arial', 'FontSize',fontSize_pt);
    set(ax.YLabel, 'FontName','Arial', 'FontSize',fontSize_pt);

    ax.XTickLabel = replaceLeadingHyphen(ax.XTickLabel);
    ax.YTickLabel = replaceLeadingHyphen(ax.YTickLabel);
end

function labels = replaceLeadingHyphen(labels)
% Replace only a leading ASCII hyphen in a tick label.
    wasCharMatrix = ischar(labels);
    if wasCharMatrix
        labels = cellstr(labels);
    elseif isstring(labels)
        labels = cellstr(labels);
    end

    minusSign = char(8722); % U+2212 mathematical minus sign
    for k = 1:numel(labels)
        if ~isempty(labels{k}) && labels{k}(1) == '-'
            labels{k}(1) = minusSign;
        end
    end

    if wasCharMatrix
        labels = char(labels);
    end
end

%% ====== Y Axis tick-control component (1-10 Natural-number coverage + Significant-digit validation) ======
function [new_y_min, new_y_max, increment] = setCustomYTicks(ax, y_min_plot, y_max_plot, ymax_data)
    y_range = y_max_plot - y_min_plot;
    if y_range <= 0 || isnan(y_range) || isinf(y_range)
        y_min_plot = -1; y_max_plot = 1; ymax_data = 0.5;
        y_range = 2;
    end
    
    % Extract the base order of magnitude
    mag = 10^floor(log10(y_range / 5)); 
    
    % Generate full coverage of 1 Up to 9 natural numbers, and intermediate values
    candidate_ds = [];
    for m = [mag/10, mag, mag*10]
        % Candidate steps covering the full range
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
        
        if first_tick > test_max
            continue;
        end
        
        ticks = first_tick : d : (test_max + d * 0.01);
        ticks = round(ticks ./ d) .* d;
        n_ticks = length(ticks);
        
        % Rule 1: Limit the number of ticks (3 ~ 8 items, Prefer 4~6)
        if n_ticks < 3 || n_ticks > 8
            continue;
        end
        
        % ★ Rule 2: Significant-digit filtering (Filter out schemes whose final digit is neither even nor 5 acceptable)
        valid_ticks_flag = true;
        for t = ticks
            v = abs(t);
            if v < 1e-10 * d, continue; end % Skip 0
            % Get significant-digit text
            s = num2str(v, '%g');
            e_idx = find(s == 'e' | s == 'E');
            if ~isempty(e_idx), s = s(1:e_idx-1); end
            s = strrep(s, '.', '');
            s = regexprep(s, '^0+', ''); % Remove leading 0
            s = regexprep(s, '0+$', ''); % Remove trailing 0
            
            % If >= 2 significant digits, Check the final digit
            if length(s) > 1
                if ~ismember(s(end), ['0','2','4','6','8','5'])
                    valid_ticks_flag = false;
                    break;
                end
            end
        end
        if ~valid_ticks_flag
            continue; % If any label is invalid for this step, reject the entire scheme
        end
        
        % Rule 3: Scoring system
        stretch = (test_max - y_max_plot) / y_range; % Stretch
        tick_penalty = max(0, n_ticks - 6) * 0.05;   % Penalty for too many labels
        
        % Limit margins to at most 25%
        top_ratio = (test_max - ymax_data) / (test_max - y_min_plot);
        top_penalty = 0;
        if top_ratio > 0.25
            top_penalty = (top_ratio - 0.25) * 50; 
        end
        
        % Prefer integer or commonly used steps
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
    
    % Fallback if all schemes are rejected
    if isempty(best_ticks)
        best_d = mag * 2; % Steps beginning with 2 are valid
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
    
    % First use MATLAB default significant-digit formatting for labels.
    % Only when nonzero ticks have unequal decimal precision, pad trailing 0; 
    % 0 Keep as "0", Retain original labels when decimal precision already matches.
    tick_labels = cell(1, length(best_ticks));
    decimal_counts = nan(1, length(best_ticks));
    for i = 1:length(best_ticks)
        if abs(best_ticks(i)) < 1e-10 * increment, best_ticks(i) = 0; end
        tick_labels{i} = num2str(best_ticks(i), '%g');

        if best_ticks(i) ~= 0
            decimal_counts(i) = localDecimalPlaces(tick_labels{i});
        end
    end

    nonzero_idx = ~isnan(decimal_counts);
    if any(nonzero_idx)
        unique_decimal_counts = unique(decimal_counts(nonzero_idx));
        if numel(unique_decimal_counts) > 1
            target_decimals = max(unique_decimal_counts);
            for i = find(nonzero_idx)
                tick_labels{i} = sprintf(['%.', num2str(target_decimals), 'f'], best_ticks(i));
            end
        end
    end
    
    set(ax, 'YTickLabel', tick_labels);
end

function n = localDecimalPlaces(label)
    % Return decimal precision in ordinary decimal notation, Support %g generated scientific notation.
    exponent_idx = find(label == 'e' | label == 'E', 1);
    if isempty(exponent_idx)
        mantissa = label;
        exponent = 0;
    else
        mantissa = label(1:exponent_idx-1);
        exponent = str2double(label(exponent_idx+1:end));
    end

    dot_idx = find(mantissa == '.', 1);
    if isempty(dot_idx)
        mantissa_decimals = 0;
    else
        mantissa_decimals = numel(mantissa) - dot_idx;
    end
    n = max(0, mantissa_decimals - exponent);
end

%% ====== Statistical utility function ======
function ranovatbl = p2Ride2Ranova(xx1)
    xx1_tmp = array2table(xx1, 'VariableNames', {'s1', 's2', 'p1', 'p2'});
    rm = fitrm(xx1_tmp, 's1,s2,p1,p2 ~ 1', ...
               'WithinDesign', table([1 -1 1 -1]', [1 1 -1 -1]', 'VariableNames', {'lambda'; 'drug'}));
    ranovatbl = ranova(rm, 'WithinModel', 'lambda*drug');
end
