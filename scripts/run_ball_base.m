%RUN_BALL_BASE Reproduce the base robust ASP instance with ball uncertainty.

repo_root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(repo_root, 'src')));
yalmip('clear');

N = 50;
data = ASPell_input(N);
param = ASPell_param(data);

[ub.schedule, ub.cost, ub.time, ub.z] = ASPell_ub_dual(data, param);
[lb.schedule, lb.cost, lb.time, lb.U] = ASPell_lb(data, param, ub);
[ex.schedule, ex.cost, ex.time, ex.wcs, ex.comp] = ...
    ASPell_ex(data, param, lb);

results.data = data;
results.upper_bound = ub;
results.lower_bound = lb;
results.exact = ex;

output_file = fullfile(repo_root, 'results', 'ball_base.mat');
% save(output_file, 'results');

fprintf('\nBase ball instance complete. Exact robust cost: %.6f\n', ex.cost);
