function data = ASPpol_input(N)
%ASPPOL_INPUT Define the budget uncertainty instance.

    nominal = 3;
    rho = 1;
    gamma = rho * N / 2.5;

    T = nominal * N;
    xbar = nominal * ones(N, 1);

    D = [ eye(N),       zeros(N);
         -eye(N),       zeros(N);
          eye(N),      -eye(N);
         -eye(N),      -eye(N);
          zeros(1, N),  ones(1, N)];

    d = [rho * ones(N, 1) + xbar;
         rho * ones(N, 1) - xbar;
         xbar;
        -xbar;
         gamma];

    data.name = "budget";
    data.N = N;
    data.D = D;
    data.d = d;
    data.T = T;
    data.co = 3;
    data.cw = 0.2;
    data.ci = 1;

end