function param = ASPell_param(data)
%ASPELL_PARAM Construct matrices for the ellipsoidal ASP formulation.

    N = data.N;
    L = size(data.Q, 2);
    M = N * (N + 5) / 2;

    %% Construct H
    blockSizes = (1:N)' + 1;
    blockSizes = [blockSizes(1:N-1); blockSizes(N-1); blockSizes(N)];
    repeatCounts = repmat(blockSizes, N + 1, 1);

    H = eye(N + 1);
    H = repelem(H(:), repeatCounts);
    H = reshape(H, M, N + 1);

    %% Construct A
    rowIdx = [0; cumsum((1:N)' + 1)] + 1;
    rowIdx(end) = rowIdx(end) - 1;

    A = zeros(M, N);

    for n = 1:N-1
        block = rot90(tril(ones(n)));
        A(rowIdx(n)+1:rowIdx(n+1)-1, 1:n) = ...
            data.cw * block;
    end

    block = tril(ones(N - 1));
    A(rowIdx(N)+1:rowIdx(N+1)-1, 1:N-1) = ...
        -data.ci * block;

    block = rot90(tril(ones(N)));
    A(rowIdx(N+1)+1:end, 1:N) = ...
        data.co * block;

    %% Cost vector
    c = [data.cw * ones(N - 1, 1);
         data.ci;
         data.co];

    %% Store parameters
    param.N = N;
    param.L = L;
    param.M = M;
    param.H = H;
    param.A = A;
    param.c = c;

end