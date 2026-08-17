function [schedule, costs, time, U_hill] = ASPell_lb(data, param, ub)
%ASPELL_LB Compute the lower bound using scenario selection and hill climbing.

    tic;

    %% Parameters
    N = data.N;
    Q = data.Q;
    xbar = data.xbar;
    rho = data.rho;

    M = param.M;

    z = ub.z;
    schedule = ub.schedule;

    sqrtQ = sqrtm(Q);

    options = sdpsettings('solver', 'mosek', 'verbose', 0);

    %% Construct initial candidate scenarios
    % The initial set contains x* and feasible ratios V_m / w_m.
    U_all = zeros(M + 1, N);
    keep = false(M + 1, 1);

    U_all(1, :) = z.x';
    keep(1) = true;

    for m = 1:M
        if z.w(m) > 1e-6
            x = z.V(m, :) / z.w(m);

            if norm(sqrtQ * (x' - xbar)) <= rho + 1e-4
                U_all(m + 1, :) = x;
                keep(m + 1) = true;
            end
        end
    end

    U_all = U_all(keep, :);

    % Remove nearly identical scenarios.
    U_all = remove_duplicates(U_all, 1e-4);

    %% Adversarial scenario selection
    C = size(U_all, 1);
    selected = false(C, 1);

    currentCost = -inf;

    while any(~selected)

        % Evaluate all unselected scenarios under the current schedule.
        scenarioCosts = -inf(C, 1);

        for i = 1:C
            if ~selected(i)
                scenarioCosts(i) = evaluate_cost( ...
                    U_all(i, :), schedule, data);
            end
        end

        % Add the scenario with the highest cost.
        [~, ind] = max(scenarioCosts);
        selected(ind) = true;

        U_selected = U_all(selected, :);
        [newSchedule, newCost] = compute_givenU( ...
            U_selected, data, options);

        % Stop when adding another scenario no longer improves the bound.
        if isfinite(currentCost) && newCost - currentCost <= 1e-3
            if nnz(selected) > 1
                selected(ind) = false;
            end
            break
        end

        currentCost = newCost;
        schedule = newSchedule;

    end

    U_adv = U_all(selected, :);
    [schedule_adv, ~] = compute_givenU(U_adv, data, options);

    %% Hill climbing
    C = size(U_adv, 1);
    U_hill = zeros(C, N);

    for i = 1:C
        U_hill(i, :) = hill_climbing( ...
            U_adv(i, :)', schedule_adv, data, param)';
    end

    %% Final lower bound
    [schedule, costs] = compute_givenU(U_hill, data, options);

    time = toc;

end


%% Hill climbing
function x_hill = hill_climbing(x, schedule, data, param)

    Q = data.Q;
    xbar = data.xbar;
    rho = data.rho;

    N = data.N;
    M = param.M;
    A = param.A;

    % From Appendix C.1: b(s) = -As.
    b = -A * schedule;

    %% Optimize w for fixed x
    idx = [0; cumsum((1:N)' + 1)] + 1;
    idx(end) = idx(end) - 1;

    w = zeros(M, 1);
    values = A * x + b;

    for i = 1:N
        [~, argmax] = max(values(idx(i) + 1:idx(i + 1)));
        w(idx(i) + argmax) = 1;
    end

    [~, argmax] = max(values(idx(N) + 1:end));
    w(idx(N) + argmax) = 1;

    %% Optimize x for fixed w
    % Closed-form solution for the ellipsoidal uncertainty set.
    direction = A' * w;

    numerator = Q \ direction;
    denominator = norm(sqrtm(Q) \ direction);

    if denominator <= 1e-12
        x_hill = x;
    else
        x_hill = xbar + rho * numerator / denominator;
    end

end


%% Solve the ASP for a finite scenario set
function [schedule, costs] = compute_givenU(U, data, options)

    N = data.N;
    T = data.T;
    L = size(U, 1);

    yalmip('clear');

    s = sdpvar(N, 1);
    tau = sdpvar(1, 1);
    w = sdpvar(L, N + 1, 'full');

    waitNext = w(:, 2:N+1);
    waitRecursion = w(:, 1:N) + U - repmat(s', L, 1);

    constraints = [
        sum(s) == T
        s >= 0
        w(:) >= 0
        w(:, 1) == 0
        waitNext(:) >= waitRecursion(:)
        data.cw * sum(w(:, 1:N), 2) ...
            + data.co * w(:, N+1) ...
            + data.ci * ( ...
                w(:, N) ...
                + sum(s(1:N-1)) ...
                - sum(U(:, 1:N-1), 2)) <= tau
    ];

    diagnostics = optimize(constraints, tau, options);

    if diagnostics.problem ~= 0
        error('Finite-scenario optimization failed: %s', ...
            diagnostics.info);
    end

    schedule = value(s);
    costs = value(tau);

end


%% Evaluate the ASP cost for one scenario
function cost = evaluate_cost(x, s, data)

    N = data.N;

    wait = zeros(N + 1, 1);

    for n = 2:N+1
        wait(n) = max( ...
            wait(n - 1) + x(n - 1) - s(n - 1), 0);
    end

    waitingCost = data.cw * sum(wait(1:N));
    overtimeCost = data.co * wait(N + 1);

    idleTime = wait(N) ...
        + sum(s(1:N-1)) ...
        - sum(x(1:N-1));

    idleCost = data.ci * idleTime;

    cost = waitingCost + idleCost + overtimeCost;

end


%% Remove nearly identical rows
function U = remove_duplicates(U, tolerance)

    keep = true(size(U, 1), 1);

    for i = 1:size(U, 1)
        if ~keep(i)
            continue
        end

        for j = i+1:size(U, 1)
            if keep(j) && norm(U(i, :) - U(j, :)) < tolerance
                keep(j) = false;
            end
        end
    end

    U = U(keep, :);

end