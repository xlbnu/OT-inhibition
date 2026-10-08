function [pl,fl]=errorCIbound(x,y,color,x_range)
mdl = fitlm(x, y); % Fit a linear model
if nargin<4
x_range = linspace(min(x), max(x), 100)';
end
[ypred, ci] = predict(mdl, x_range); % Predictions and confidence intervals

% Plot scatter points, Regression line and shaded confidence interval
% Fill confidence interval(Shading between two regression lines)
hold on
fl=fill([x_range; flipud(x_range)], [ci(:,1); flipud(ci(:,2))], ...
    color, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', '95% 置信区间');
% Plot regression line
pl=plot(x_range, ypred, 'color',color, 'LineWidth', 1, 'DisplayName', '回归线');%
end