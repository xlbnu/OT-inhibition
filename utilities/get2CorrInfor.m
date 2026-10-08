function get2CorrInfor(x, y)

result = corr_bootstrap(x, y, ...
    'Type', 'Pearson', ...
    'NBoot', 100000, ...
    'Alpha', 0.05, ...
    'Seed', 42);

% p < 0.0001 Use scientific notation when, Retain 2 significant digits
if result.p_boot < 0.001
    pText = sprintf('%.2e', result.p_boot);
else
    pText = sprintf('%.4f', result.p_boot);
end

fprintf('\nr = %.3g, p = %s, 95%% CI = [%.3g, %.3g]\n\n', ...
    result.r, pText, result.ci(1), result.ci(2));

end