function [schedule, costs, time, z] = ASPpol_ub_dual(data, param)
%ASPPOL_UB_DUAL Compute the RPT upper bound using the dual formulation.

    tic;

    %% Parameters
    N = data.N;
    D = data.D;
    d = data.d;
    T = data.T;

    N2 = param.N2;
    M = param.M;
    H = param.H;
    A = param.A;

    %% Decision variables
    yalmip('clear');

    y = sdpvar(1, 1);
    w = sdpvar(M, 1);
    x = sdpvar(N2, 1);
    V = sdpvar(M, N2);

    %% Constraints
    theta = -A' * w;

    constraints = [
        theta + y * ones(N2, 1) >= 0
        D * x <= d
        w >= 0
        H' * w == ones(N + 1, 1)
        w * d' - V * D' >= 0
        H' * V == repmat(x', N + 1, 1)
    ];

    %% Objective and optimization
    objective = T * y - trace(A' * V);

    options = sdpsettings( ...
        'solver', 'gurobi', ...
        'verbose', 0);

    diagnostics = optimize(constraints, objective, options);

    if diagnostics.problem ~= 0
        error('Upper-bound optimization failed: %s', ...
            diagnostics.info);
    end

    %% Store output
    dualValues = dual(constraints(1));

    schedule = dualValues(1:N);
    costs = -value(objective);

    z.V = value(V);
    z.w = value(w);
    z.x = value(x);

    time = toc;

end