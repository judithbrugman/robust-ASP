%RUN_BUDGET_SCALABILITY Compute bounds and runtimes for varying N.

repo_root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(repo_root, 'src')));
yalmip('clear');

N_values = 2:100;
output_file = fullfile(fileparts(mfilename('fullpath')), '..', 'results', ...
    'budget_scale.mat');

results = struct([]);

for k = 1:numel(N_values)
    N = N_values(k);
    fprintf('\n=== N = %d ===\n', N);

    data = ASPpol_input(N);
    param = ASPpol_param(data);

    [ub.schedule, ub.cost, ub.time, ub.z] = ASPpol_ub_dual(data, param);
    % [ub_support.schedule, ub_support.cost, ub_support.time, ub_support.z] = ...
    %     ASPpol_ub_support(data, param);
    [lb.schedule, lb.cost, lb.time, lb.U] = ASPpol_lb(data, param, ub);
    [ex.schedule, ex.cost, ex.time, ex.wcs, ex.comp] = ...
        ASPpol_ex(data, param, lb);

    results(k).N = N;
    results(k).upper_bound = ub;
    % results(k).upper_bound_support = ub_support;
    results(k).lower_bound = lb;
    results(k).exact = ex;

    % save(output_file, 'results');

    fprintf('Finished N = %d: UB = %.4f, LB = %.4f, exact = %.4f\n', ...
        N, ub.cost, lb.cost, ex.cost);

end
