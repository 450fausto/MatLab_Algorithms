function respuesta = sistema_33(w)
% n-butanol + Agua, NRTL
% x* = 0.017489391493040133
% f* = -0.021368188827368576
%
% Malinen, I., Kangas, J., & Tanskanen, J. (2012).
% A new Newton homotopy based method for the robust determination of all the stationary points of the tangent plane distance function.
% Chemical Engineering Science, 84, 266–275. doi:10.1016/j.ces.2012.08.037

    % ---- Parámetros del problema ----
    cantidad_variables = 2; %#ok<NASGU>  % fidelidad al Fortran
    z = [0.3, 0.7];

    % Parámetros NRTL (tau y G)
    tau = [ ...
        0.0,       0.00745177; ...
        3.80209,   0.0        ];
    G = [ ...
        1.0,       0.99702373; ...
        0.21852912,1.0        ];

    % ---- Construcción de x a partir de w ----
    % Dos componentes: x = [x1, x2] con x2 = 1 - x1.
    x = [ w(1), 1 - w(1) ];

    % ---- ln(gamma) con NRTL ----
    lngx = ln_gamma_nrtl(x, tau, G);
    lngz = ln_gamma_nrtl(z, tau, G);

    % Alinear dimensiones por si ln_gamma_nrtl devuelve columna
    if size(lngx,1) > 1, lngx = lngx.'; end
    if size(lngz,1) > 1, lngz = lngz.'; end

    % ---- Penalización por fracciones negativas (emula PACK de Fortran) ----
    negativos = x(x < 0);
    criterio  = sum(negativos);   % <= 0; <0 si existe algún negativo

    % ---- Función objetivo de estabilidad (tangente de Gibbs) ----
    if criterio < 0
        respuesta = -criterio * 1.0e6;   % 10.0E5 en Fortran
    else
        % Nota: fidelidad al Fortran; x == 0 no se penaliza explícitamente.
        respuesta = sum( x .* ( log(x) + lngx - log(z) - lngz ) );
    end
end
