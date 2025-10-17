function respuesta = sistema_08(w)
% etilenglicol (1), lauril alcohol (2) y nitrometano (3).
% Etilenglicol + dodecanol + nitrometano, UNIQUAC
% x* = 0.75425313547529549, 0.0022192732704503513
% f* = -0.11395145263903501
%
% Tessier S. R., Brennecke J. F. and Stadtherr M. A. (2000).
% Reliable phase stability analysis for excess Gibbs energy models.
% Chemical Engineering Science, 55(10), 1785-1796. Doi: 10.1016/s0009-2509(99)00442-x
% Problem 4
%
% McDonald C. M. and Floudas C. A. (1994).
% Global Optimization for the phase stability problem.
% AIChE Journal, 41(7), 1798-1814. Doi: 10.1002/aic.690410715
% Example 5

    % ---- Parámetros del problema ----
    cantidad_variables = 3;

    % Composición global
    z = [0.4, 0.3, 0.3];

    % Parámetros UNIQUAC
    r  = [2.4088, 8.8495, 2.0086];
    q  = [2.248,  7.372,  1.868 ];
    qq = [2.248,  7.372,  1.868 ];

    % Matriz tau (tal como en el Fortran)
    tau = [ ...
        1.0,     0.432589, 0.830749; ...
        0.789593, 1.0,     0.354992; ...
        0.204736, 0.636678, 1.0     ];

    % ---- Construcción de x a partir de w ----
    % Tres componentes: x = [x1, x2, x3] con x3 = 1 - x1 - x2.
    x = [ w(1), w(2), 1 - w(1) - w(2) ];

    % ---- ln(gamma) con UNIQUAC ----
    lngx = ln_gamma_uniquac(x, r, q, qq, tau);
    lngz = ln_gamma_uniquac(z, r, q, qq, tau);

    % Alinear dimensiones a fila por si ln_gamma_uniquac devuelve columna
    if size(lngx,1) > 1, lngx = lngx.'; end
    if size(lngz,1) > 1, lngz = lngz.'; end

    % ---- Penalización por fracciones negativas ----
    % Emulamos PACK de Fortran: sumar sólo los elementos negativos de x.
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
