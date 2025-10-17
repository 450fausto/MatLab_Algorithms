function respuesta = sistema_07(w)
% Tolueno + agua + anilina, NRTL
% x* = 0.000066937752536364210, 0.99686525905361334
% f* = -0.29453555322324276
%
% McDonald C. M. and Floudas C. A. (1994).
% Global Optimization for the phase stability problem.
% AIChE Journal, 41(7), 1798-1814. Doi: 10.1002/aic.690410715
% Example 2

    % ---- Parámetros del problema ----
    cantidad_variables = 3;

    % Composición global z
    z = [0.29989, 0.20006, 0.50005];

    % Parámetros NRTL (tau y G)
    tau = [ ...
        0.0,     4.93035, 1.59806; ...
        7.77063, 0.0,     4.18462; ...
        0.03509, 1.27932, 0.0     ];
    G = [ ...
        1.0,    0.29370, 0.61914; ...
        0.14500, 1.0,    0.23984; ...
        0.98953, 0.64629, 1.0    ];

    % ---- Construcción de x a partir de w ----
    % Tres componentes: x = [x1, x2, x3] con x3 = 1 - x1 - x2.
    x = [ w(1), w(2), 1 - w(1) - w(2) ];

    % ---- ln(gamma) con NRTL ----
    lngx = ln_gamma_nrtl(x, tau, G);
    lngz = ln_gamma_nrtl(z, tau, G);

    % Alinear dimensiones a fila por si ln_gamma_nrtl devuelve columna
    if size(lngx,1) > 1, lngx = lngx.'; end
    if size(lngz,1) > 1, lngz = lngz.'; end

    % ---- Penalización por fracciones negativas ----
    % Emulamos PACK de Fortran: sumar sólo los elementos negativos de x.
    negativos = x(x < 0);
    criterio  = sum(negativos);   % <= 0; <0 si existe algún negativo

    % ---- Función objetivo de estabilidad (tangente de Gibbs) ----
    if criterio < 0
        respuesta = -criterio * 1.0e6;   % 10.0E5
    else
        % Nota: se mantiene fidelidad al Fortran; x == 0 no se penaliza explícitamente.
        respuesta = sum( x .* ( log(x) + lngx - log(z) - lngz ) );
    end
end
