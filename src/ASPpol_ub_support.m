function [schedule, costs, time, z] = ASPpol_ub_support(data, param)
%ASPPOL_UB_SUPPORT Compute the RPT upper bound using the support function.

    tic;

    %% Parameters
    N = data.N;
    D = data.D;
    d = data.d;
    T = data.T;

    N2 = param.N2;
    K = param.K;
    M = param.M;
    H = param.H;
    A = param.A;

    %% Decision variables
    yalmip('clear');

    w1 = sdpvar(K, 1);
    w2 = sdpvar(M, 1);
    w3 = sdpvar(N + 1, 1);
    W4 = sdpvar(M, K, 'full');
    W5 = sdpvar(N + 1, N2, 'full');

    s = sdpvar(N, 1);
    tau = sdpvar(1, 1);

    % Extend the schedule with zeros when the polyhedral representation
    % contains auxiliary uncertainty variables.
    sExt = [s; zeros(N2 - N, 1)];

    %% Construct b(s)
    b = -A * sExt;

    %% Constraints
    c1 = d' * w1 + ones(1, N + 1) * w3 <= tau;
    c2 = -w2 + H * w3 - W4 * d == b;
    c3 = D' * w1 - W5' * ones(N + 1, 1) == zeros(N2, 1);
    c4 = W4 * D + H * W5 == A;

    constraints = [
        c1
        c2
        c3
        c4
        w1 >= 0
        w2 >= 0
        W4(:) >= 0
        sum(s) == T
        s >= 0
    ];

    %% Optimize
    options = sdpsettings( ...
        'solver', 'mosek', ...
        'verbose', 0);

    diagnostics = optimize(constraints, tau, options);

    if diagnostics.problem ~= 0
        error('Support-function upper-bound optimization failed: %s', ...
            diagnostics.info);
    end

    %% Store output
    schedule = value(s);
    costs = value(tau);

    z.w = -dual(c2);
    z.x = -dual(c3);
    z.V = -dual(c4);

    time = toc;

end