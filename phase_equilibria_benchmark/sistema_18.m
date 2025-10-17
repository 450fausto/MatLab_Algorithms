function respuesta = sistema_18(w)
% Acetonitrilo + benceno + nC7, NRTL
% x* = 0.91143780627404147, 0.023584205816360025
% f* = -0.10845798046210155
%
% Gecegormez, H., & Demirel, Y. (2005).
% Phase stability analysis using interval Newton method with NRTL model.
% Fluid Phase Equilibria, 237(1-2), 48–58. doi:10.1016/j.fluid.2005.08.014
% Problem 11

    % ---- Parámetros del problema ----
    cantidad_variables = 3;

    % Composición global (z)
    z = [0.4, 0.05, 0.55];

    % Parámetros NRTL (tau y G)
    tau = [ ...
        0.0,        0.5661821, 2.3187177; ...
        0.4472257,  0.0,       1.4918912; ...
        0.6964173, -0.5982783, 0.0       ];
    G = [ ...
        1.0,       0.9175650, 0.7589032; ...
        0.9343016, 1.0,       0.6290071; ...
        0.9204803, 1.2043230, 1.0       ];

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
        respuesta = -criterio * 1.0e6;   % 10.0E5 en Fortran
    else
        % Nota: se mantiene fidelidad al Fortran; x == 0 no se penaliza explícitamente.
        respuesta = sum( x .* ( log(x) + lngx - log(z) - lngz ) );
    end
end
