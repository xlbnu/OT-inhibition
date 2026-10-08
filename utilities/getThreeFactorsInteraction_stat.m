
function getThreeFactorsInteraction_stat(xx1)
w = [-1 1 1 -1 1 -1 -1 1];

interaction_score = xx1 * w';

[h, p, ci, stats] = ttest(interaction_score);

fprintf('三因素交互：t(%d) = %.3f, p = %.4f\n', ...
        stats.df, stats.tstat, p);

% Equivalent repeated-measures ANOVA results
F = stats.tstat^2;
fprintf('等价表示：F(1,%d) = %.3f, p = %.4f\n', ...
        stats.df, F, p);

% Within-participant standardized effect size Cohen''s dz
dz = mean(interaction_score) / std(interaction_score);
% fprintf('Cohen''s dz = %.3f\n', dz);
%
% data: n×8, One participant per row, One experimental condition per column

% Convert the matrix into a table
varNames = {'Y111','Y112','Y121','Y122', ...
            'Y211','Y212','Y221','Y222'};

T = array2table(xx1, 'VariableNames', varNames);

% Define 8 within-participant factor levels for each column
% Row order must match data column order exactly
A = categorical([1;1;1;1;2;2;2;2]);% A B
context = categorical([1;1;2;2;1;1;2;2]);% social nonsocial
drug = categorical([1;2;1;2;1;2;1;2]);% OT PL

withinDesign = table(A, context, drug);

disp(withinDesign);

% Fit the repeated-measures model
% No between-participant factors, Therefore the between-participant model contains only an intercept: ~ 1
rm = fitrm(T, 'Y111-Y222 ~ 1', ...
           'WithinDesign', withinDesign);

% Three-factor repeated-measures ANOVA
anovaTable = ranova(rm, 'WithinModel', 'A*context*drug');

disp(anovaTable);


% Three-way interaction row
rowNames = string(anovaTable.Properties.RowNames);
idxABC = contains(rowNames, 'A:context:drug') & ...
         ~contains(rowNames, 'Error');

% Three-way interaction error row
idxErrorABC = contains(rowNames, 'A:context:drug') & ...
              contains(rowNames, 'Error');
F_ABC = anovaTable.F(idxABC);
df1 = anovaTable.DF(idxABC);
df2 = anovaTable.DF(idxErrorABC);

eta_p2_F = (F_ABC * df1) / (F_ABC * df1 + df2);

fprintf('\nA×context×drug：F(%g,%g)=%.3f, p=%f, partical η²=%.3f (Cohen''s d = %.3f)\n', ...
    df1, df2, F_ABC, anovaTable.pValue(15), eta_p2_F,dz);
end
