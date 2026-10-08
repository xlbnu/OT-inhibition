function RemoveBoxplotlines(gca,lineWith,lineColor)
if nargin<2
lineWith=1.5;
end
if nargin<3
lineColor='k';
end
% Remove end caps on upper and lower whiskers
lav = findobj(gca, 'Tag', 'Lower Adjacent Value'); % Get the lower end-cap handle
uav = findobj(gca, 'Tag', 'Upper Adjacent Value'); % Get the upper end-cap handle
set([lav, uav], 'LineStyle', 'none'); % Hide both
% Make upper and lower whiskers solid
% By default, boxplot whiskers may be dashed
lw = findobj(gca, 'Tag', 'Lower Whisker'); % Find lower whiskers
uw = findobj(gca, 'Tag', 'Upper Whisker'); % Find upper whiskers
set([lw, uw], 'LineStyle', '-'); % Set solid line style

% Thicken median lines
median_line = findobj(gca, 'Tag', 'Median'); % Get median-line handles
set(median_line, 'LineWidth', lineWith); % Set line width(Larger values give thicker lines, Usually defaults to 1)
% Optional: Adjust additional properties(Such as color)
set(median_line, 'Color', lineColor); % Set median lines to black
end