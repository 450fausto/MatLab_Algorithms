function respuesta = sistema_39(w)
% SBA + DSBE + water, UNIQUAC
% C3H8O + C4H10O + H2O
% x* = 0.034449071209442503, 0.000000049173133921076047
% f* = -0.073510203500150875
%
% McDonald C. M. and Floudas C. A. (1994).
% Global Optimization for the phase stability problem.
% AIChE Journal, 41(7), 1798-1814. Doi: 10.1002/aic.690410715
% Example 6

    % ---- Parámetros del problema ----
    cantidad_variables = 3; %#ok<NASGU>

    % Composición global (z)
    z = [0.40, 0.04, 0.56];

    % Parámetros UNIQUAC
    r  = [3.9235, 6.0909, 0.9200];
    q  = [3.6640, 5.1680, 1.4000];
    qq = [4.0643, 5.7409, 1.6741];

    % Matriz tau (tal como en el Fortran)
    tau = [ ...
        1.0,        1.31226645, 0.55571771; ...
        0.5620479,  1.0,        0.64605486; ...
        0.8660343,  0.00436255, 1.0       ];

    % ---- Construcción de x a partir de w ----
    % Tres componentes: x = [x1, x2, x3] con x3 = 1 - x1 - x2.
    x = [ w(1), w(2), 1 - w(1) - w(2) ];

    % ---- ln(gamma) con UNIQUAC ----
    lngx = ln_gamma_uniquac(x, r, q, qq, tau);
    lngz = ln_gamma_uniquac(z, r, q, qq, tau);

    % Alinear dimensiones a fila por si ln_gamma_uniquac devuelve columna
    if size(lngx,1) > 1, lngx = lngx.'; end
    if size(lngz,1) > 1, lngz = lngz.'; end

    % ---- Penalización por fracciones negativas (emula PACK de Fortran) ----
    negativos = x(x < 0);
    criterio  = sum(negativos);   % < 0 si existe algún negativo

    % ---- Función objetivo de estabilidad (tangente de Gibbs) ----
    if criterio < 0
        respuesta = -criterio * 1.0e6;   % 10.0E5 en Fortran
    else
        % Nota: fidelidad al Fortran; x == 0 no se penaliza explícitamente.
        respuesta = sum( x .* ( log(x) + lngx - log(z) - lngz ) );
    end
end
