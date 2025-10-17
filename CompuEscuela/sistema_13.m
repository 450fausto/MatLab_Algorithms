function respuesta = sistema_13(w)
% n-Propanol + n-butanol + benceno + agua, NRTL
% x* = 0.018114780903023944, 0.00061999525113539625, 0.0044844004034368520
% f* = -0.33982213819701940
%
% Tessier S. R., Brennecke J. F. and Stadtherr M. A. (2000).
% Reliable phase stability analysis for excess Gibbs energy models.
% Chemical Engineering Science, 55(10), 1785-1796. Doi: 10.1016/s0009-2509(99)00442-x
% Problem 2

    % ---- Parámetros del problema ----
    cantidad_variables = 4;

    % Composición global (z)
    z = [0.148, 0.052, 0.60, 0.20];

    % Parámetros NRTL (tau y G)
    tau = [ ...
        0.0,     2.16486,  0.23689,  0.13060; ...
       -1.2007,  0.0,     -0.09730,  0.19154; ...
        2.01911, 1.73912,  0.0,      4.01932; ...
        2.31985, 4.31706,  4.09334,  0.0     ];
    G = [ ...
        1.0,     0.34320,  0.93449,  0.96384; ...
        1.80967, 1.0,      1.02932,  0.93623; ...
        0.56132, 0.59659,  1.0,      0.32322; ...
        0.51986, 0.22649,  0.31656,  1.0     ];

    % ---- Construcción de x a partir de w ----
    % Cuatro componentes: x = [x1, x2, x3, x4] con x4 = 1 - x1 - x2 - x3.
    x = [ w(1), w(2), w(3), 1 - w(1) - w(2) - w(3) ];

    % ---- ln(gamma) con NRTL ----
    lngx = ln_gamma_nrtl(x, tau, G);
    lngz = ln_gamma_nrtl(z, tau, G);

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
