function respuesta = sistema_26(w)
% Hexanol + nitrometano + agua, UNIQUAC
% x* = 0.00059220030402009878, 0.025320498362062040
% f* = -0.17182033237996530
%
% Negahban, S., Willhite, G. P., & Wales, S. M. (1988).
% Modeling of Three-Phase Liquid/Liquid Equilibria.
% SPE Reservoir Engineering, 3(03), 1017–1024. doi:10.2118/14936-pa

    % ---- Parámetros del problema ----
    cantidad_variables = 3;

    % Composición global (z)
    z = [0.2, 0.4, 0.4];

    % Parámetros UNIQUAC
    r  = [4.8031, 2.0086, 0.92 ];
    q  = [4.132,  1.868,  1.4  ];
    qq = [4.132,  1.868,  1.4  ];

    % Matriz tau (tal como en el Fortran)
    tau = [ ...
        1.0,        0.34788273, 0.77732613; ...
        0.85240166, 1.0,        0.24767237; ...
        0.33360179, 0.50876062, 1.0       ];

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
    % Emula PACK de Fortran: sumar solo los elementos negativos de x.
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
