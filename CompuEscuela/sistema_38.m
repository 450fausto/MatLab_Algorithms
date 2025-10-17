function respuesta = sistema_38(w)
% Toluene + water, UNIQUAC
% x* = 0.00049200239761908724
% f* = -0.61887272125625625
%
% McDonald C. M. and Floudas C. A. (1994).
% Global Optimization for the phase stability problem.
% AIChE Journal, 41(7), 1798-1814. Doi: 10.1002/aic.690410715
% Example 4

    % ---- Parámetros del problema ----
    cantidad_variables = 2; %#ok<NASGU>
    z = [0.5, 0.5];

    % Parámetros UNIQUAC
    r  = [3.92, 0.92];
    q  = [2.97, 1.40];
    qq = [2.97, 1.00];   % Asimetría superficial efectiva (como en el Fortran)

    % Matriz tau (tal como en el Fortran)
    tau = [ ...
        1.0,    0.09867; ...
        0.59673, 1.0    ];

    % ---- Construcción de x a partir de w ----
    % Dos componentes: x = [x1, x2] con x2 = 1 - x1.
    x = [ w(1), 1 - w(1) ];

    % ---- ln(gamma) con UNIQUAC ----
    lngx = ln_gamma_uniquac(x, r, q, qq, tau);
    lngz = ln_gamma_uniquac(z, r, q, qq, tau);

    % Alinear dimensiones por si ln_gamma_uniquac devuelve columna
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
