function respuesta = sistema_25(w)
% Agua + 1,4-dicianobutano, NRTL
% x* = 0.99603654451877344
% f* = -0.14491953271850552
%
% Gecegormez, H., & Demirel, Y. (2005).
% Phase stability analysis using interval Newton method with NRTL model.
% Fluid Phase Equilibria, 237(1-2), 48�58. doi:10.1016/j.fluid.2005.08.014
% Problem 8

    % ---- Par�metros del problema ----
    cantidad_variables = 2;
    z = [0.4, 0.6];

    % Par�metros NRTL (tau y G) tal como en el Fortran
    tau = [ ...
        0.0,      4.6707725; ...
        0.6918617, 0.0      ];
    G = [ ...
        1.0,      0.2462934; ...
        0.8125657, 1.0      ];

    % ---- Construcci�n de x a partir de w ----
    % Dos componentes: x = [x1, x2] con x2 = 1 - x1.
    x = [ w(1), 1 - w(1) ];

    % ---- ln(gamma) con NRTL ----
    lngx = ln_gamma_nrtl(x, tau, G);
    lngz = ln_gamma_nrtl(z,  tau, G);

    % Alinear dimensiones por si ln_gamma_nrtl devuelve columna
    if size(lngx,1) > 1, lngx = lngx.'; end
    if size(lngz,1) > 1, lngz = lngz.'; end

    % ---- Penalizaci�n por fracciones negativas (emula PACK de Fortran) ----
    negativos = x(x < 0);
    criterio  = sum(negativos);   % <= 0; <0 si existe alg�n negativo

    % ---- Funci�n objetivo de estabilidad (tangente de Gibbs) ----
    if criterio < 0
        respuesta = -criterio * 1.0e6;   % 10.0E5 en Fortran
    else
        % Nota: fidelidad al Fortran; x == 0 no se penaliza expl�citamente.
        respuesta = sum( x .* ( log(x) + lngx - log(z) - lngz ) );
    end
end
