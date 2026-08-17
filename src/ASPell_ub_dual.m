function [schedule, costs, time, z] = ASPell_ub_dual(data, param)
%ASPELL_UB_DUAL Compute the RPT upper bound using the dual formulation.

    tic;

    %% Parameters
    N = data.N;
    Q = data.Q;
    xbar = data.xbar;
    rho = data.rho;
    T = data.T;

    M = param.M;
    H = param.H;
    A = param.A;

    sqrtQ = sqrtm(Q);

    %% Decision variables
    yalmip('clear');

    y = sdpvar(1, 1);
    w = sdpvar(M, 1);
    x = sdpvar(N, 1);
    V = sdpvar(M, N);

    %% Constraints
    theta = -A' * w;

    constraints = [
        theta + y * ones(N, 1) >= 0
        norm(sqrtQ * (x - xbar)) <= rho
        w >= 0
        H' * w == ones(N + 1, 1)
        H' * V == repmat(x', N + 1, 1)
    ];

    for m = 1:M
        constraints = [
            constraints
            norm(sqrtQ * (V(m, :)' - w(m) * xbar)) <= rho * w(m)
        ];
    end

    %% Objective and optimization
    objective = T * y - trace(A' * V);

    options = sdpsettings( ...
        'solver', 'mosek', ...
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