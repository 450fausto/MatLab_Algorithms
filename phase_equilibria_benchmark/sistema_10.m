function respuesta = sistema_10(w)
% n-Propanol + n-butanol + agua, NRTL
% x* = 0.0094074726576420242, 0.019052055523897715
% f* = -0.011609328088991122
%
% Tessier S. R., Brennecke J. F. and Stadtherr M. A. (2000).
% Reliable phase stability analysis for excess Gibbs energy models.
% Chemical Engineering Science, 55(10), 1785-1796. Doi: 10.1016/s0009-2509(99)00442-x
% Problem 1
%
% McDonald C. M. and Floudas C. A. (1994).
% Global Optimization for the phase stability problem.
% AIChE Journal, 41(7), 1798-1814. Doi: 10.1002/aic.690410715
% Example 3

    % ---- Parámetros del problema ----
    cantidad_variables = 3;

    % Composición global z
    z = [0.04, 0.16, 0.80];

    % Parámetros NRTL (tau y G) tal como en el Fortran
    tau = [ ...
        0.0,     -0.61259, -0.07149; ...
        0.71640,  0.0,      0.90047; ...
        2.7425,   3.51307,  0.0     ];
    G = [ ...
        1.0,       1.2017478, 1.021678; ...
        0.8066060, 1.0,       0.6490629; ...
        0.4392221, 0.1852084, 1.0      ];

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
