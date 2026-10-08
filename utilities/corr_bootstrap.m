function result = corr_bootstrap(x, y, varargin)
%CORR_BOOTSTRAP Correlation with a paired nonparametric bootstrap test.
%
%   RESULT = CORR_BOOTSTRAP(X,Y) calculates Pearson's correlation and uses
%   10,000 paired bootstrap samples to obtain a 95% percentile confidence
%   interval and an approximate two-sided test of H0: rho = 0.
%
%   RESULT = CORR_BOOTSTRAP(...,'Type','Spearman','NBoot',5000,
%       'Alpha',0.05,'Seed',1,'Plot',true)
%
%   X and Y must be real vectors of equal length. Observations for which
%   either value is NaN are removed as pairs. Each bootstrap sample resamples
%   the remaining (x,y) pairs, thereby preserving their joint structure.
%
%   Name-value options:
%       Type  - 'Pearson' (default) or 'Spearman'
%       NBoot - number of bootstrap samples (default 10000)
%       Alpha - significance level (default 0.05)
%       Seed  - scalar random seed; [] leaves RNG unchanged (default [])
%       Plot  - show bootstrap distribution (default false)
%
%   RESULT fields include r, ci, p_boot, reject, n, n_boot_valid and r_boot.
%   The p-value is based on the bootstrap distribution centered at zero.

    p = inputParser;
    p.FunctionName = mfilename;
    addRequired(p, 'x', @(z) isvector(z) && isnumeric(z) && isreal(z));
    addRequired(p, 'y', @(z) isvector(z) && isnumeric(z) && isreal(z));
    addParameter(p, 'Type', 'Pearson', @(s) ischar(s) || (isstring(s) && isscalar(s)));
    addParameter(p, 'NBoot', 10000, @(z) isscalar(z) && z >= 100 && z == fix(z));
    addParameter(p, 'Alpha', 0.05, @(z) isscalar(z) && z > 0 && z < 1);
    addParameter(p, 'Seed', [], @(z) isempty(z) || (isscalar(z) && isnumeric(z) && isfinite(z)));
    addParameter(p, 'Plot', false, @(z) isscalar(z) && (islogical(z) || isnumeric(z)));
    parse(p, x, y, varargin{:});

    type = validatestring(p.Results.Type, {'Pearson','Spearman'});
    x = double(x(:));
    y = double(y(:));
    if numel(x) ~= numel(y)
        error('corr_bootstrap:SizeMismatch', 'X and Y must have the same number of elements.');
    end

    keep = ~isnan(x) & ~isnan(y);
    x = x(keep);
    y = y(keep);
    n = numel(x);
    if n < 4
        error('corr_bootstrap:TooFewPairs', 'At least four complete observation pairs are required.');
    end
    if any(~isfinite(x)) || any(~isfinite(y))
        error('corr_bootstrap:NonFinite', 'X and Y may contain NaN, but not Inf or -Inf.');
    end

    rObserved = localCorrelation(x, y, type);
    if isnan(rObserved)
        error('corr_bootstrap:ConstantInput', 'Correlation is undefined because an input is constant.');
    end

    oldRng = [];
    if ~isempty(p.Results.Seed)
        oldRng = rng;
        rng(p.Results.Seed, 'twister');
    end
    cleanup = onCleanup(@() restoreRng(oldRng)); %#ok<NASGU>

    nBoot = p.Results.NBoot;
    rBoot = nan(nBoot, 1);
    for b = 1:nBoot
        idx = randi(n, n, 1);
        rBoot(b) = localCorrelation(x(idx), y(idx), type);
    end
    rBoot = rBoot(isfinite(rBoot));
    if numel(rBoot) < max(100, 0.9*nBoot)
        error('corr_bootstrap:TooFewValidSamples', ...
            'Too many bootstrap samples had undefined correlations.');
    end

    alpha = p.Results.Alpha;
    ci = localPercentile(rBoot, [alpha/2, 1-alpha/2]);

    % Center the empirical bootstrap distribution to approximate H0: rho=0.
    rNull = rBoot - mean(rBoot);
    pBoot = (sum(abs(rNull) >= abs(rObserved)) + 1) / (numel(rNull) + 1);

    result = struct( ...
        'r', rObserved, ...
        'ci', ci, ...
        'p_boot', pBoot, ...
        'reject', pBoot < alpha, ...
        'alpha', alpha, ...
        'type', type, ...
        'n', n, ...
        'n_removed', numel(keep) - n, ...
        'n_boot_requested', nBoot, ...
        'n_boot_valid', numel(rBoot), ...
        'r_boot', rBoot);

    if logical(p.Results.Plot)
        figure('Color', 'w');
        histogram(rBoot, 'Normalization', 'pdf', 'FaceColor', [0.25 0.55 0.85]);
        hold on;
        yl = ylim;
        plot([rObserved rObserved], yl, 'r-', 'LineWidth', 2);
        plot([ci(1) ci(1)], yl, 'k--', 'LineWidth', 1.5);
        plot([ci(2) ci(2)], yl, 'k--', 'LineWidth', 1.5);
        plot([0 0], yl, ':', 'Color', [0.25 0.25 0.25], 'LineWidth', 1.2);
        hold off;
        xlabel(sprintf('%s correlation', type));
        ylabel('Bootstrap density');
        title(sprintf('r = %.3f, %.1f%% CI [%.3f, %.3f], p = %.4g', ...
            rObserved, 100*(1-alpha), ci(1), ci(2), pBoot));
        legend('Bootstrap estimates', 'Observed r', 'CI limits', '', 'Zero', ...
            'Location', 'best');
        box off;
    end
end

function r = localCorrelation(x, y, type)
    if strcmpi(type, 'Spearman')
        x = localTiedRank(x);
        y = localTiedRank(y);
    end
    x = x - mean(x);
    y = y - mean(y);
    denom = sqrt(sum(x.^2) * sum(y.^2));
    if denom == 0
        r = NaN;
    else
        r = sum(x .* y) / denom;
        r = max(-1, min(1, r));
    end
end

function ranks = localTiedRank(values)
    [sorted, order] = sort(values);
    ranks = zeros(size(values));
    first = 1;
    n = numel(values);
    while first <= n
        last = first;
        while last < n && sorted(last + 1) == sorted(first)
            last = last + 1;
        end
        ranks(order(first:last)) = (first + last) / 2;
        first = last + 1;
    end
end

function q = localPercentile(values, probabilities)
    values = sort(values(:));
    n = numel(values);
    q = zeros(size(probabilities));
    for k = 1:numel(probabilities)
        pos = 1 + (n - 1) * probabilities(k);
        lo = floor(pos);
        hi = ceil(pos);
        q(k) = values(lo) + (pos - lo) * (values(hi) - values(lo));
    end
end

function restoreRng(oldRng)
    if ~isempty(oldRng)
        rng(oldRng);
    end
end
