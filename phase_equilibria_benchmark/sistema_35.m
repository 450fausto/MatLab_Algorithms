function respuesta = sistema_35(w)
% metil ter-butil éter (MTBE) + metanol + Agua, NRTL
% x* = 0.0097753532660367921, 0.12602139563983467
% f* = -0.091365474230144267
%
% Arce, A., Blanco, A., Blanco, M., Soto, A. and Vidal, I. (1994),
% Liquid-liquid equilibria of water + methanol + (MTBE or TAME) mixtures.
% Can. J. Chem. Eng., 72: 935-938. https://doi.org/10.1002/cjce.5450720522

    % ---- Parámetros del problema ----
    cantidad_variables = 3; %#ok<NASGU>  % fidelidad al Fortran

    % Composición global (z)
    z = [0.39726, 0.20292, 0.399824];

    % Parámetros NRTL (tau y G)
    tau = [ ...
         0.0,        -1.17786349,  1.21200738; ...
         1.42391414,  0.0,         -2.66362569; ...
         5.6528593,    3.58141875,  0.0        ];
    G = [ ...
        1.0,          1.26563339,  0.78474106; ...
        0.75217759,   1.0,         1.70356845; ...
        0.32284858,   0.48856451,  1.0        ];

    % ---- Construcción de x a partir de w ----
    % Tres componentes: x = [x1, x2, x3] con x3 = 1 - x1 - x2.
    x = [ w(1), w(2), 1 - w(1) - w(2) ];

    % ---- ln(gamma) con NRTL ----
    lngx = ln_gamma_nrtl(x, tau, G);
    lngz = ln_gamma_nrtl(z, tau, G);

    % Alinear dimensiones a fila por si ln_gamma_nrtl devuelve columna
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
