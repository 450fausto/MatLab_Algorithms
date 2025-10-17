function respuesta = sistema_16(y)
% C1 + CO2 + nC6 + H2S, SRK
% x* = 0.96094306587817768, 0.013921929386259638, 0.0000099597629233678039
% f* = -0.066114216390560962
%
% Ivanov, B. B., Galushko, A. A., & Stateva, R. P. (2013).
% Phase Stability Analysis with Equations of State—A Fresh Look from  a Different Perspective.
% Industrial & Engineering Chemistry Research, 52(32), 11208–11223. doi:10.1021/ie401072x
% Problem 7

    % ---- Parámetros del problema ----
    cantidad_variables = 4;
    P = 42.5;
    T = 200.0;

    % Composición global (z)
    z = [0.5, 0.0574, 0.0263, 0.4163];

    % Propiedades críticas y factor acéntrico
    Tc = [190.6, 304.2, 507.5, 373.2];
    Pc = [46.00, 73.80, 30.10, 89.40];
    w  = [0.008, 0.225, 0.299, 0.100];

    % Matriz de interacción binaria k_ij
    k = [ ...
        0.0,    0.0933, 0.0360, 0.0800; ...
        0.0933, 0.0,    0.1180, 0.0989; ...
        0.0360, 0.1180, 0.0,    0.0500; ...
        0.0800, 0.0989, 0.0500, 0.0    ];

    % ---- Construcción de x a partir de y ----
    % 3 variables libres: y(1:3); la cuarta se cierra por suma = 1.
    x = zeros(1, cantidad_variables);
    x(1:cantidad_variables-1) = [y(1), y(2), y(3)];
    x(cantidad_variables)     = 1 - sum(x);

    % ---- ln(phi) con SRK ----
    ln_phi_x = ln_phi_srk(x, Tc, Pc, w, k, P, T);
    ln_phi_z = ln_phi_srk(z, Tc, Pc, w, k, P, T);

    % Alinear dimensiones por si ln_phi_srk devuelve columna
    if size(ln_phi_x,1) > 1, ln_phi_x = ln_phi_x.'; end
    if size(ln_phi_z,1) > 1, ln_phi_z = ln_phi_z.'; end

    % ---- Penalización por fracciones negativas ----
    negativos = x(x < 0);
    criterio  = sum(negativos);   % <= 0; <0 si existe algún negativo

    % ---- Función objetivo de estabilidad (tangente de Gibbs) ----
    if criterio < 0
        respuesta = -criterio * 1.0e6;   % 10.0E5 en Fortran
    else
        % Nota: fidelidad al Fortran; x == 0 no se penaliza explícitamente.
        respuesta = sum( x .* ( log(x) + ln_phi_x - log(z) - ln_phi_z ) );
    end
end
