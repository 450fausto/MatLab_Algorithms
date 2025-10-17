function respuesta = sistema_04(w)
% n-Acetato de butilo + agua NRTL
% x* = 0.0042086615362195347 
% f* = -0.032638313157465078
%
% McDonald C. M. and Floudas C. A. (1994).
% Global Optimization for the phase stability problem.
% AIChE Journal, 41(7), 1798-1814. Doi: 10.1002/aic.690410715
% Example 1

    % ---- Parámetros ----
    cantidad_variables = 2;
    z = [0.5, 0.5];

    % Parámetros NRTL (tau y G) como en el Fortran
    tau = [ ...
        0.0,     3.00498; ...
        4.69071, 0.0     ];
    G = [ ...
        1.0,     0.30800; ...
        0.15909, 1.0     ];

    % ---- Construcción de x a partir de w ----
    x = [w(1), 1 - w(1)];

    % ---- ln(gamma) con NRTL ----
    lngx = ln_gamma_nrtl(x, tau, G);
    lngz = ln_gamma_nrtl(z, tau, G);

    % Alinear dimensiones por si ln_gamma_nrtl devuelve columna
    if size(lngx,1) > 1, lngx = lngx.'; end
    if size(lngz,1) > 1, lngz = lngz.'; end

    % ---- Penalización por fracciones negativas ----
    negativos = x(x < 0);
    criterio  = sum(negativos);   % <= 0

    % ---- Función objetivo (tangente de Gibbs) ----
    if criterio < 0
        respuesta = -criterio * 1.0e6;   % 10.0E5
    else
        % Nota: se mantiene el comportamiento de Fortran (no se penaliza x == 0)
        respuesta = sum( x .* ( log(x) + lngx - log(z) - lngz ) );
    end
end
