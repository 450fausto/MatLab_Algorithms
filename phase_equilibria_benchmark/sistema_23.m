function respuesta = sistema_23(w)
% Agua + ácido cítrico + 2-butanol, NRTL
% x* = 0.083523360596904961, 0.0070909368991194941
% f* = -0.021738486434602477
%
% Gecegormez, H., & Demirel, Y. (2005).
% Phase stability analysis using interval Newton method with NRTL model.
% Fluid Phase Equilibria, 237(1-2), 48–58. doi:10.1016/j.fluid.2005.08.014
% Problem 12

    % ---- Parámetros del problema ----
    cantidad_variables = 3;

    % Composición global (z)
    z = [0.15, 0.10, 0.75];

    % Parámetros NRTL (tau y G)
    tau = [ ...
         0.0,         0.9889317,  2.9732685; ...
        13.7521382,   0.0,       -1.3581754; ...
         0.5249036,   7.4341774,  0.0       ];
    G = [ ...
        1.0,          0.6887706,  0.2472300; ...
        0.0056008823, 1.0,        1.3199010; ...
        0.7813714,    0.2188763,  1.0       ];

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
