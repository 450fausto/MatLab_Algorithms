function respuesta = sistema_05(y)
% C1 + C2 + C3 + iC4 + C4 + iC5 + C5 + C6 + iC15, SRK
% x* = 0.93301293603872271, 0.030185805930259917, 0.018572154791578286, 0.0054796462078695657, 0.0074918054834172730, 0.00061840983302726218, 0.0038403394484845519, 0.00079713025666419972
% x* = 0.946076,            0.043568,             0.007851,             0.000675,              0.001247,              0.000197,               0.000264,              0.000121, 0.000000
% f* = -1.4571297453173282
%
% Rangaiah, G. P. (2001).
% Evaluation of genetic algorithms and simulated annealing for phase equilibrium and stability problems.
% Fluid Phase Equilibria, 187-188, 83–109. doi:10.1016/s0378-3812(01)00528-3
%
% Ivanov, B. B., Galushko, A. A., & Stateva, R. P. (2013).
% Phase Stability Analysis with Equations of State—A Fresh Look from  a Different Perspective.
% Industrial & Engineering Chemistry Research, 52(32), 11208–11223. doi:10.1021/ie401072x
% Problem 10

    % ---- Parámetros del problema ----
    cantidad_variables = 9;
    P = 20.1;
    T = 314;

    % Composición global (z)
    z  = [0.614, 0.10259, 0.04985, 0.008989, 0.02116, 0.00722, 0.01187, 0.01435, 0.16998];
    %            C1      C2       C3       iC4       C4       iC5      C5       C6       iC15

    % Propiedades críticas y factor acéntrico
    Tc = [190.6, 305.4, 369.8, 408.1, 425.2, 460.4, 469.6, 507.4, 707.0];
    Pc = [46.0,  48.84, 42.46, 36.48, 38.0,  33.84, 33.74, 29.69, 15.2];
    w  = [0.008, 0.098, 0.152, 0.176, 0.193, 0.227, 0.251, 0.296, 0.706];

    % Matriz de interacción binaria k_ij (aquí toda en cero)
    k = zeros(cantidad_variables);

    % ---- Construcción de x a partir de y ----
    % 8 variables libres: y(1:8); la novena se cierra por suma = 1.
    x = zeros(1, cantidad_variables);
    x(1:cantidad_variables-1) = [y(1), y(2), y(3), y(4), y(5), y(6), y(7), y(8)];
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
        respuesta = -criterio * 1.0e6;     % 10.0E5 en Fortran
    else
        % Nota: se conserva la fidelidad al Fortran; x == 0 no se penaliza explícitamente.
        respuesta = sum( x .* ( log(x) + ln_phi_x - log(z) - ln_phi_z ) );
    end
end
