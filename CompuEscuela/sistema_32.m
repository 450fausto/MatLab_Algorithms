function respuesta = sistema_32(w)
% Agua + ácido cítrico, NRTL
% x* = 0.4
% f* = 0.0
%
% Gecegormez, H., & Demirel, Y. (2005).
% Phase stability analysis using interval Newton method with NRTL model.
% Fluid Phase Equilibria, 237(1-2), 48–58. doi:10.1016/j.fluid.2005.08.014
% Problem 5

    % ---- Parámetros del problema ----
    cantidad_variables = 2; %#ok<NASGU>  % fidelidad al Fortran
    z = [0.4, 0.6];

    tau = [ ...
        0.0,        0.9889317; ...
        13.7521382, 0.0       ];
    G = [ ...
        1.0,         0.6887706; ...
        0.0056008823,1.0      ];

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
