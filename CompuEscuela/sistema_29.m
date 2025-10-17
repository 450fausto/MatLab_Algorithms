function respuesta = sistema_29(w)
% n-Pentanol + 2,2-dimetilbutano, NRTL
% x* = 0.15003426724818841
% f* = -0.000098192104616469118
%
% Gecegormez, H., & Demirel, Y. (2005).
% Phase stability analysis using interval Newton method with NRTL model.
% Fluid Phase Equilibria, 237(1-2), 48–58. doi:10.1016/j.fluid.2005.08.014
% Problem 1

    % ---- Parámetros del problema ----
    cantidad_variables = 2; %#ok<NASGU>  % mantenido por fidelidad al Fortran
    z = [0.1, 0.9];

    % Parámetros NRTL (tau y G)
    tau = [ ...
        0.0,       1.0114044; ...
        2.9413053, 0.0       ];
    G = [ ...
        1.0,       0.5527024; ...
        0.1782906, 1.0       ];

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
