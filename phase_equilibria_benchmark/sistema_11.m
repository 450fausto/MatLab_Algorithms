function respuesta = sistema_11(y)
% C1 + CO2 + H2S, SRK
% x* = 0.90389409900234807, 0.038708064233775250
% f* = -0.19727432314258303
%
% Stateva, R. P., & Tsvetkov, S. G. (1994).
% A diverse approach for the solution of the isothermal multiphase flash problem. Application to vapor-liquid-liquid systems.
% The Canadian Journal of Chemical Engineering, 72(4), 722–734. doi:10.1002/cjce.5450720422
%
% Ivanov, B. B., Galushko, A. A., & Stateva, R. P. (2013).
% Phase Stability Analysis with Equations of State—A Fresh Look from  a Different Perspective.
% Industrial & Engineering Chemistry Research, 52(32), 11208–11223. doi:10.1021/ie401072x
% Problem 6

    % ---- Parámetros del problema ----
    cantidad_variables = 3;
    P = 48.60;
    T = 227.55;

    % Composición global (z)
    z = [0.4989, 0.0988, 0.4023];

    % Propiedades críticas y factor acéntrico
    Tc = [190.564, 304.1282, 373.1];
    Pc = [45.992,   73.773,   90.0 ];
    w  = [0.01142,  0.22394,  0.1005];

    % Matriz k_ij
    k = [ ...
        0.0,    0.0933, 0.0823; ...
        0.0933, 0.0,    0.0989; ...
        0.0823, 0.0989, 0.0    ];

    % ---- Construcción de x a partir de y ----
    % Tres componentes: x = [x1, x2, x3] con x3 = 1 - x1 - x2.
    x = [ y(1), y(2), 1 - y(1) - y(2) ];

    % ---- ln(phi) con SRK ----
    ln_phi_x = ln_phi_srk(x, Tc, Pc, w, k, P, T);
    ln_phi_z = ln_phi_srk(z, Tc, Pc, w, k, P, T);

    % Alinear dimensiones por si ln_phi_srk devuelve columna
    if size(ln_phi_x,1) > 1, ln_phi_x = ln_phi_x.'; end
    if size(ln_phi_z,1) > 1, ln_phi_z = ln_phi_z.'; end

    % ---- Penalización por fracciones negativas ----
    % Emulamos PACK de Fortran: sumar sólo los elementos negativos de x.
    negativos = x(x < 0);
    criterio  = sum(negativos);   % <= 0; <0 si existe algún negativo

    % ---- Función objetivo de estabilidad (tangente de Gibbs) ----
    if criterio < 0
        respuesta = -criterio * 1.0e6;   % 10.0E5
    else
        % Nota: se mantiene fidelidad al Fortran; x == 0 no se penaliza explícitamente.
        respuesta = sum( x .* ( log(x) + ln_phi_x - log(z) - ln_phi_z ) );
    end
end
