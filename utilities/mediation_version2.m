
% [results] = bootstrap_mediation_detailed((x_gaba), (x_lambda), (x_v), 100000);

% function [results] = bootstrap_mediation_detailed(A, B, C, n_iterations)
function results = mediation_version2(A, B, C, n_iterations)
% MEDIATION_VERSION2 Bootstrap Mediation analysis (A -> B -> C)
%
% Input:
%   A: Independent variable (Nx1 vector)
%   B: Mediator (Nx1 vector)
%   C: Dependent variable (Nx1 vector)
%   n_iterations: Bootstrap Number of iterations(Recommended 10000)
%
% Output:
%   results: Contains path coefficients, standard errors, t value, degrees of freedom, Bootstrap
%            confidence intervals and p values in a structure
%
% Note:
%   a, b, c', c Path t values and degrees of freedom come from OLS regression.
%   ab Indirect-effect t values use Sobel/Delta method approximation; 
%   Primary inference for the indirect effect remains based on Bootstrap confidence intervals.

%% 0. Organize and validate inputs

A = A(:);
B = B(:);
C = C(:);

if ~(length(A) == length(B) && length(B) == length(C))
    error('A、B、C 的长度必须相同。');
end

if ~isscalar(n_iterations) || n_iterations < 1 || ...
        n_iterations ~= floor(n_iterations)
    error('n_iterations 必须是正整数。');
end

% Remove any observations containing NaN or Inf
valid_idx = isfinite(A) & isfinite(B) & isfinite(C);
A = A(valid_idx);
B = B(valid_idx);
C = C(valid_idx);

n = length(A);

if n <= 3
    error('删除缺失值后，至少需要 4 个有效观测。');
end

if std(A) == 0 || std(B) == 0 || std(C) == 0
    error('A、B、C 均必须具有非零方差。');
end

% Standardize
A = zscore(A);
B = zscore(B);
C = zscore(C);

%% 1. Original-sample path coefficients, standard errors, t values and degrees of freedom

% ---------- path a: A -> B ----------
X_a = [ones(n, 1), A];
beta_a = X_a \ B;

path_a = beta_a(2);

resid_a = B - X_a * beta_a;
df_a = n - rank(X_a);
mse_a = sum(resid_a.^2) / df_a;
cov_beta_a = mse_a * pinv(X_a' * X_a);

se_a = sqrt(cov_beta_a(2, 2));
t_a = path_a / se_a;
p_t_a = 2 * tcdf(-abs(t_a), df_a);

% ---------- path b and c': A, B -> C ----------
X_bc = [ones(n, 1), A, B];
beta_bc = X_bc \ C;

path_c_prime = beta_bc(2);
path_b = beta_bc(3);

resid_bc = C - X_bc * beta_bc;
df_bc = n - rank(X_bc);
mse_bc = sum(resid_bc.^2) / df_bc;
cov_beta_bc = mse_bc * pinv(X_bc' * X_bc);

se_c_prime = sqrt(cov_beta_bc(2, 2));
se_b = sqrt(cov_beta_bc(3, 3));

t_c_prime = path_c_prime / se_c_prime;
t_b = path_b / se_b;

df_c_prime = df_bc;
df_b = df_bc;

p_t_c_prime = 2 * tcdf(-abs(t_c_prime), df_c_prime);
p_t_b = 2 * tcdf(-abs(t_b), df_b);

% ---------- path c: A -> C ----------
X_c = [ones(n, 1), A];
beta_c = X_c \ C;

path_c = beta_c(2);

resid_c = C - X_c * beta_c;
df_c = n - rank(X_c);
mse_c = sum(resid_c.^2) / df_c;
cov_beta_c = mse_c * pinv(X_c' * X_c);

se_c = sqrt(cov_beta_c(2, 2));
t_c = path_c / se_c;
p_t_c = 2 * tcdf(-abs(t_c), df_c);

% ---------- Indirect effect ab ----------
path_ab = path_a * path_b;

% Sobel/Delta method standard errors
se_ab = sqrt(path_b^2 * se_a^2 + path_a^2 * se_b^2);
t_ab = path_ab / se_ab;

% ab Not a single OLS regression coefficient, No strictly unique residual degrees of freedom.
% Use conservative approximate degrees of freedom here.
df_ab = min(df_a, df_b);
p_t_ab = 2 * tcdf(-abs(t_ab), df_ab);

%% 2. Bootstrap

boot_a = zeros(n_iterations, 1);
boot_b = zeros(n_iterations, 1);
boot_c_prime = zeros(n_iterations, 1);
boot_c = zeros(n_iterations, 1);
boot_ab = zeros(n_iterations, 1);

for i = 1:n_iterations

    idx = randi(n, n, 1);

    A_b = A(idx);
    B_b = B(idx);
    C_b = C(idx);

    % path a
    beta_a_b = [ones(n, 1), A_b] \ B_b;
    boot_a(i) = beta_a_b(2);

    % path b and c'
    beta_bc_b = [ones(n, 1), A_b, B_b] \ C_b;
    boot_c_prime(i) = beta_bc_b(2);
    boot_b(i) = beta_bc_b(3);

    % path c
    beta_c_b = [ones(n, 1), A_b] \ C_b;
    boot_c(i) = beta_c_b(2);

    % Indirect effect
    boot_ab(i) = boot_a(i) * boot_b(i);
end

%% 3. Bootstrap 95% Percentile confidence intervals

ci_a = prctile(boot_a, [2.5, 97.5]);
ci_b = prctile(boot_b, [2.5, 97.5]);
ci_c_prime = prctile(boot_c_prime, [2.5, 97.5]);
ci_c = prctile(boot_c, [2.5, 97.5]);
ci_ab = prctile(boot_ab, [2.5, 97.5]);

%% 4. Bootstrap Two-tailed p value

% Add 1 correction, Avoid exactly p = 0
calc_pval = @(dist) min(1, ...
    2 * (min(sum(dist <= 0), sum(dist >= 0)) + 1) / ...
    (length(dist) + 1));

p_boot_a = calc_pval(boot_a);
p_boot_b = calc_pval(boot_b);
p_boot_c_prime = calc_pval(boot_c_prime);
p_boot_c = calc_pval(boot_c);
p_boot_ab = calc_pval(boot_ab);

%% 5. Print results

% fprintf('\n');
% fprintf('================================================================================================================\n');
% fprintf('                     Bootstrap Mediation analysis results: A -> B -> C(N = %d, Number of iterations = %d)\n', ...
%     n, n_iterations);
% fprintf('================================================================================================================\n');
% fprintf('%-12s %10s %10s %10s %8s %22s %12s %12s\n', ...
%     'path', 'Effect estimate', 'standard errors', 't value', 'df', ...
%     'Bootstrap 95% CI', 't test p', 'Bootstrap p');
% fprintf('----------------------------------------------------------------------------------------------------------------\n');
% 
% fprintf('%-12s %10.4f %10.4f %10.4f %8d (%9.4f, %9.4f) %12.4g %12.4g\n', ...
%     'Path a', path_a, se_a, t_a, df_a, ...
%     ci_a(1), ci_a(2), p_t_a, p_boot_a);
% 
% fprintf('%-12s %10.4f %10.4f %10.4f %8d (%9.4f, %9.4f) %12.4g %12.4g\n', ...
%     'Path b', path_b, se_b, t_b, df_b, ...
%     ci_b(1), ci_b(2), p_t_b, p_boot_b);
% 
% fprintf('%-12s %10.4f %10.4f %10.4f %8d (%9.4f, %9.4f) %12.4g %12.4g\n', ...
%     'Path c''', path_c_prime, se_c_prime, t_c_prime, df_c_prime, ...
%     ci_c_prime(1), ci_c_prime(2), p_t_c_prime, p_boot_c_prime);
% 
% fprintf('%-12s %10.4f %10.4f %10.4f %8d (%9.4f, %9.4f) %12.4g %12.4g\n', ...
%     'Path c', path_c, se_c, t_c, df_c, ...
%     ci_c(1), ci_c(2), p_t_c, p_boot_c);
% 
% fprintf('%-12s %10.4f %10.4f %10.4f %8d (%9.4f, %9.4f) %12.4g %12.4g\n', ...
%     'Path ab', path_ab, se_ab, t_ab, df_ab, ...
%     ci_ab(1), ci_ab(2), p_t_ab, p_boot_ab);
% 
% fprintf('================================================================================================================\n');
% fprintf('Note: Path ab 's t values and degrees of freedom are Sobel/Delta method approximate.\n');
% fprintf('    For indirect-effect inference, prioritize Bootstrap 95%% CI and Bootstrap p values.\n');
% 
% if ci_ab(1) > 0 || ci_ab(2) < 0
%     fprintf('Conclusion: Indirect effect is significant, Bootstrap 95%% CI Excludes 0.\n');
% else
%     fprintf('Conclusion: Indirect effect is not significant, Bootstrap 95%% CI including 0.\n');
% end
% 
% fprintf('================================================================================================================\n\n');

%% 5. Print results

fprintf('\n============================================================\n');
fprintf('Bootstrap 中介分析结果（N = %d，迭代次数 = %d）\n', ...
    n, n_iterations);
fprintf('============================================================\n');

fprintf('Path a（自变量→中介）: β=%.4f, t(%d)=%.4f, 95%%CI=[%.4f, %.4f], p=%.4g\n', ...
    path_a, df_a, t_a, ci_a(1), ci_a(2), p_boot_a);

fprintf('Path b（中介→因变量）: β=%.4f, t(%d)=%.4f, 95%%CI=[%.4f, %.4f], p=%.4g\n', ...
    path_b, df_b, t_b, ci_b(1), ci_b(2), p_boot_b);

fprintf('Path c''（直接效应）: β=%.4f, t(%d)=%.4f, 95%%CI=[%.4f, %.4f], p=%.4g\n', ...
    path_c_prime, df_c_prime, t_c_prime, ...
    ci_c_prime(1), ci_c_prime(2), p_boot_c_prime);

fprintf('Path c（总效应）: β=%.4f, t(%d)=%.4f, 95%%CI=[%.4f, %.4f], p=%.4g\n', ...
    path_c, df_c, t_c, ci_c(1), ci_c(2), p_boot_c);

fprintf('Path ab（间接效应）: β=%.4f, t(NA)=NA, 95%%CI=[%.4f, %.4f], p=%.4g\n', ...
    path_ab, ci_ab(1), ci_ab(2), p_boot_ab);

fprintf('============================================================\n');
fprintf('注：95%%CI 为 Bootstrap 百分位置信区间；p 为 Bootstrap p 值。\n');
fprintf('    ab 为路径系数乘积，不适用常规回归 t 值和自由度。\n\n');
%% 6. Return results

results.n = n;
results.n_iterations = n_iterations;

results.a.estimate = path_a;
results.a.se = se_a;
results.a.t = t_a;
results.a.df = df_a;
results.a.p_t = p_t_a;
results.a.ci_boot = ci_a;
results.a.p_boot = p_boot_a;

results.b.estimate = path_b;
results.b.se = se_b;
results.b.t = t_b;
results.b.df = df_b;
results.b.p_t = p_t_b;
results.b.ci_boot = ci_b;
results.b.p_boot = p_boot_b;

results.c_prime.estimate = path_c_prime;
results.c_prime.se = se_c_prime;
results.c_prime.t = t_c_prime;
results.c_prime.df = df_c_prime;
results.c_prime.p_t = p_t_c_prime;
results.c_prime.ci_boot = ci_c_prime;
results.c_prime.p_boot = p_boot_c_prime;

results.c.estimate = path_c;
results.c.se = se_c;
results.c.t = t_c;
results.c.df = df_c;
results.c.p_t = p_t_c;
results.c.ci_boot = ci_c;
results.c.p_boot = p_boot_c;

results.ab.estimate = path_ab;
results.ab.se_sobel = se_ab;
results.ab.t_sobel = t_ab;
results.ab.df_approx = df_ab;
results.ab.p_sobel = p_t_ab;
results.ab.ci_boot = ci_ab;
results.ab.p_boot = p_boot_ab;

% Optional: Saving Bootstrap distributions, For subsequent plotting or inspection
results.bootstrap.a = boot_a;
results.bootstrap.b = boot_b;
results.bootstrap.c_prime = boot_c_prime;
results.bootstrap.c = boot_c;
results.bootstrap.ab = boot_ab;

end