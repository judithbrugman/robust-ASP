function param = ASPpol_param(data)
%ASPPOL_PARAM Construct model parameters for the polyhedral ASP.

    N = data.N;
    N2 = size(data.D, 2);
    K = size(data.D, 1);
    M = N * (N + 5) / 2;

    %% Construct H
    H = eye(N + 1);

    blockSizes = (1:N)' + 1;
    blockSizes = [blockSizes(1:N-1); blockSizes(N-1); blockSizes(N)];

    repetitions = repmat(blockSizes, N + 1, 1);
    H = repelem(H(:), repetitions);
    H = reshape(H, [M, N + 1]);

    %% Construct A and c
    idx = [0; cumsum((1:N)' + 1)] + 1;
    idx(end) = idx(end) - 1;

    A = zeros(M, N2);
    c = [repmat(data.cw, N - 1, 1); data.ci; data.co];

    % Waiting-time terms
    for i = 1:N-1
        mat = rot90(tril(ones(i)));
        A(idx(i) + 1:idx(i + 1) - 1, 1:i) = ...
            A(idx(i) + 1:idx(i + 1) - 1, 1:i) + data.cw * mat;
    end

    % Idle-time terms
    mat = tril(ones(N - 1));
    A(idx(N) + 1:idx(N + 1) - 1, 1:N-1) = -data.ci * mat;

    % Overtime terms
    mat = rot90(tril(ones(N)));
    A(idx(N + 1) + 1:end, 1:N) = data.co * mat;

    %% Store parameters
    param.N2 = N2;
    param.K = K;
    param.M = M;
    param.H = H;
    param.c = c;
    param.A = A;

end