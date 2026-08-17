function [schedule, costs, time, U_wc, comp_wc] = ASPell_ex(data, param, lb)
%ASPELL_EX Compute the exact robust solution using the cutting-set method.

    tic;

    %% Parameters
    tolerance = 1e-2;
    options = sdpsettings('solver', 'mosek', 'verbose', 0);

    %% Initialize with lower-bound scenario set
    U = lb.U;
    schedule = lb.schedule;
    lowerBound = lb.cost;

    %% Cutting-set algorithm
    while true

        % Find worst-case scenario for the current schedule.
        [xWorst, upperBound] = solve_worst_case( ...
            schedule, data, param, options);

        gap = upperBound - lowerBound;

        if gap <= tolerance
            break
        end

        % Add worst-case scenario.
        U = [U; xWorst'];

        % Re-optimize the schedule over the enlarged scenario set.
        [schedule, lowerBound] = compute_givenU(U, data, options);

    end

    %% Store results
    costs = lowerBound;
    [U_wc, comp_wc] = find_worstcase(schedule, costs, U, data);

    time = toc;

end


%% Worst-case scenario for a fixed schedule
function [xWorst, costWorst] = solve_worst_case( ...
        schedule, data, param, options)

    N = param.N;
    Q = data.Q;
    xbar = data.xbar;
    rho = data.rho;

    sqrtQ = sqrtm(Q);

    % Construct a valid big-M from the ellipsoidal uncertainty set.
    Qinv = Q \ eye(N);
    radius = rho * sqrt(max(diag(Qinv), 0));

    x_upper = xbar + radius;
    x_lower = xbar - radius;

    bigM = max([ ...
        sum(max(x_upper, 0)), ...
        data.T - min(x_lower), ...
        data.T]);

    yalmip('clear');

    x = sdpvar(N, 1);
    wait = sdpvar(N + 1, 1);
    z = binvar(N, 1);

    objective = ...
        data.cw * sum(wait(1:N)) ...
        + data.ci * ( ...
            wait(N) ...
            + sum(schedule(1:N-1)) ...
            - sum(x(1:N-1))) ...
        + data.co * wait(N+1);

    constraints = [
        norm(sqrtQ * (x - xbar)) <= rho
        wait(1) == 0
        wait(2:N+1) >= 0
        wait(2:N+1) >= wait(1:N) + x - schedule
        wait(2:N+1) <= bigM * (1 - z)
        wait(2:N+1) <= wait(1:N) + x - schedule + bigM * z
    ];

    diagnostics = optimize(constraints, -objective, options);

    if diagnostics.problem ~= 0
        error('Worst-case optimization failed: %s', ...
            diagnostics.info);
    end

    xWorst = value(x);
    costWorst = value(objective);

end


%% Solve the ASP for a finite scenario set
function [schedule, costs] = compute_givenU(U, data, options)

    N = data.N;
    L = size(U, 1);

    yalmip('clear');

    s = sdpvar(N, 1);
    tau = sdpvar(1, 1);
    wait = sdpvar(L, N + 1, 'full');

    waitNext = wait(:, 2:N+1);
    waitRecursion = wait(:, 1:N) + U - repmat(s', L, 1);

    constraints = [
        sum(s) == data.T
        s >= 0
        wait(:) >= 0
        wait(:, 1) == 0
        waitNext(:) >= waitRecursion(:)
        data.cw * sum(wait(:, 1:N), 2) ...
            + data.ci * ( ...
                wait(:, N) ...
                + sum(s(1:N-1)) ...
                - sum(U(:, 1:N-1), 2)) ...
            + data.co * wait(:, N+1) <= tau
    ];

    diagnostics = optimize(constraints, tau, options);

    if diagnostics.problem ~= 0
        error('Finite-scenario optimization failed: %s', ...
            diagnostics.info);
    end

    schedule = value(s);
    costs = value(tau);

end


%% Identify worst-case scenarios in the final scenario set
function [U_wc, comp] = find_worstcase( ...
        schedule, robustCost, U, data)

    N = data.N;
    C = size(U, 1);

    isWorstCase = false(C, 1);
    comp = zeros(C, 3);

    for i = 1:C

        x = U(i, :)';
        wait = zeros(N + 1, 1);

        for n = 2:N+1
            wait(n) = max( ...
                wait(n-1) + x(n-1) - schedule(n-1), ...
                0);
        end

        waitingCost = data.cw * sum(wait(1:N));

        idleCost = data.ci * ( ...
            wait(N) ...
            + sum(schedule(1:N-1)) ...
            - sum(x(1:N-1)));

        overtimeCost = data.co * wait(N+1);

        scenarioCost = waitingCost + idleCost + overtimeCost;

        comp(i, :) = [waitingCost, idleCost, overtimeCost];

        if scenarioCost >= robustCost - 1e-3
            isWorstCase(i) = true;
        end

    end

    U_wc = U(isWorstCase, :);
    comp = comp(isWorstCase, :);

end