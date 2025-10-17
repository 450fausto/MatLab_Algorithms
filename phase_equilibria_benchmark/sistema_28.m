function respuesta = sistema_28(w)
% ter-amil-metil éter (TAME) + metanol + Agua, UNIQUAC
% x* = 0.017316622873774425, 0.19943335247020999
% f* = -0.054947505251734662
%
% Arce, A., Blanco, A., Blanco, M., Soto, A. and Vidal, I. (1994),
% Liquid-liquid equilibria of water + methanol + (MTBE or TAME) mixtures.
% Can. J. Chem. Eng., 72: 935-938. https://doi.org/10.1002/cjce.5450720522

    % ---- Parámetros del problema ----
    cantidad_variables = 3;

    % Composición global (z)
    z = [0.2933, 0.2685, 0.4382];

    % Parámetros UNIQUAC
    r  = [4.7422, 1.4311, 0.92 ];
    q  = [4.172,  1.432,  1.4  ];
    qq = [4.172,  1.432,  1.4  ];

    % Matriz tau (tal como en el Fortran)
    tau = [ ...
        1.0,        0.88199016, 0.11659849; ...
        1.39626166, 1.0,        2.60702299; ...
        0.76820694, 1.18991148, 1.0       ];
    % Alternativa con más decimales (comentada en el Fortran):
    % tau = [ ...
    %     1.0,              0.8819901636917687, 0.11659849060794848; ...
    %     1.396261661823313, 1.0,               2.607022994551858;  ...
    %     0.7682069422235138,1.1899114826234083,1.0 ];

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
    criterio  = sum(negativos);   % <= 0; <0 si existe algún negativo

    % ---- Función objetivo de estabilidad (tangente de Gibbs) ----
    if criterio < 0
        respuesta = -criterio * 1.0e6;   % 10.0E5 en Fortran
    else
        % Nota: fidelidad al Fortran; x == 0 no se penaliza explícitamente.
        respuesta = sum( x .* ( log(x) + lngx - log(z) - lngz ) );
    end
end
