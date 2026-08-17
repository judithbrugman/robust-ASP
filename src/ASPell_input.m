function data = ASPell_input(N)
%ASPELL_INPUT Define the ball uncertainty instance.

    L = 10.63;
    nominal = 3;

    xbar = nominal * ones(N, 1);
    Q = eye(N);
    rho = sqrt(L);
    T = nominal * N;

    data.name = "ball";
    data.N = N;
    data.Q = Q;
    data.xbar = xbar;
    data.rho = rho;
    data.T = T;
    data.co = 3;
    data.cw = 0.2;
    data.ci = 1;

end