%RUN_BUDGET_BASE Reproduce the base robust ASP instance with budget uncertainty.

repo_root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(repo_root, 'src')));
yalmip('clear');

N = 30;
data = ASPpol_input(N);
param = ASPpol_param(data);

[ub.schedule, ub.cost, ub.time, ub.z] = ASPpol_ub_dual(data, param);
[ub_support.schedule, ub_support.cost, ub_support.time, ub_support.z] = ...
        ASPpol_ub_support(data, param);
[lb.schedule, lb.cost, lb.time, lb.U] = ASPpol_lb(data, param, ub);
[ex.schedule, ex.cost, ex.time, ex.wcs, ex.comp] = ...
    ASPpol_ex(data, param, lb);

results.data = data;
results.upper_bound = ub;
results.lower_bound = lb;
results.exact = ex;

output_file = fullfile(repo_root, 'results', 'budget_base.mat');
% save(output_file, 'results');

fprintf('\nBase budget instance complete. Exact robust cost: %.6f\n', ex.cost);
