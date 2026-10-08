function setAxesForPPT(ax, fontSize_pt, axesLineWidth_pt)
% Apply final axis styling and use U+2212 for negative tick labels.
% Positive tick labels are left unchanged.
set(ax, ...
    'FontName','Arial', ...
    'FontSize',fontSize_pt, ...
    'LineWidth',axesLineWidth_pt, ...
    'TickLabelInterpreter','none');

set(ax.XLabel, ...
    'FontName','Arial', ...
    'FontSize',fontSize_pt, ...
    'Interpreter','none');

set(ax.YLabel, ...
    'FontName','Arial', ...
    'FontSize',fontSize_pt, ...
    'Interpreter','none');

if ~isempty(ax.Title)
    set(ax.Title, ...
        'FontName','Arial', ...
        'FontSize',fontSize_pt, ...
        'Interpreter','none');
end
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
