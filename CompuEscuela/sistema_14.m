function respuesta = sistema_14(w)
% n-Propanol + n-butanol + benceno + etanol + agua, NRTL
% x* = 0.024314644551625304, 0.00054503387520804633, 0.0017265770398075587, 0.035514954118328224, 0.57410251482629460
% f* = -0.10430190637039251
%
% Tessier S. R., Brennecke J. F. and Stadtherr M. A. (2000).
% Reliable phase stability analysis for excess Gibbs energy models.
% Chemical Engineering Science, 55(10), 1785-1796. Doi: 10.1016/s0009-2509(99)00442-x
% Problem 3

    % ---- Parámetros del problema ----
    cantidad_variables = 5;

    % Composición global (z)
    z = [0.148, 0.052, 0.50, 0.10, 0.20];

    % Parámetros NRTL (tau y G) tal como en el Fortran
    tau = [ ...
        0.0,      2.16486,  0.23686,  3.78001,  0.13060; ...
       -1.20070,  0.0,     -0.09730, -1.15187, -0.20374; ...
        2.01911,  1.73912,  0.0,      1.85228,  3.73758; ...
       -0.10979,  1.16315,  0.47676,  0.0,     -0.14651; ...
        2.31985,  5.22337,  6.45226,  2.17820,  0.0     ];
    G = [ ...
        1.0,      0.34320,  0.93450,  0.35902,  0.96384; ...
        3.91873,  1.0,      1.02931,  1.41931,  1.06238; ...
        0.56132,  0.59670,  1.0,      0.57907,  0.36864; ...
        1.03030,  0.70216,  0.86880,  1.0,      1.04035; ...
        0.51986,  0.21196,  0.17857,  0.55537,  1.0     ];

    % ---- Construcción de x a partir de w ----
    % Cinco componentes: x = [x1, x2, x3, x4] y x5 = 1 - x1 - x2 - x3 - x4.
    x = [ w(1), w(2), w(3), w(4), 1 - w(1) - w(2) - w(3) - w(4) ];

    % ---- ln(gamma) con NRTL ----
    lngx = ln_gamma_nrtl(x, tau, G);
    lngz = ln_gamma_nrtl(z,  tau, G);

    % Alinear dimensiones a fila por si ln_gamma_nrtl devuelve columna
    if size(lngx,1) > 1, lngx = lngx.'; end
    if size(lngz,1) > 1, lngz = lngz.'; end

    % ---- Penalización por fracciones negativas ----
    negativos = x(x < 0);
    criterio  = sum(negativos);   % <= 0; <0 si existe algún negativo

    % ---- Función objetivo de estabilidad (tangente de Gibbs) ----
    if criterio < 0
        respuesta = -criterio * 1.0e6;   % 10.0E5 en Fortran
    else
        % Nota: se mantiene fidelidad al Fortran; x == 0 no se penaliza explícitamente.
        respuesta = sum( x .* ( log(x) + lngx - log(z) - lngz ) );
    end
end
