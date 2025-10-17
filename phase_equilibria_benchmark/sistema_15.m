function respuesta = sistema_15(y)
% C1 + C2 + C3 + nC4 + nC5 + N2, SRK
% x* = 0.51369609732804555, 0.036763277102425343, 0.012830038702293621, 0.0050251743678840196, 0.0033672222544560114
% f* = -0.0015129260816979101
%
% Stateva, R. P., & Tsvetkov, S. G. (1994).
% A diverse approach for the solution of the isothermal multiphase flash problem. Application to vapor-liquid-liquid systems.
% The Canadian Journal of Chemical Engineering, 72(4), 722–734. doi:10.1002/cjce.5450720422
%
% Ivanov, B. B., Galushko, A. A., & Stateva, R. P. (2013).
% Phase Stability Analysis with Equations of State—A Fresh Look from  a Different Perspective.
% Industrial & Engineering Chemistry Research, 52(32), 11208–11223. doi:10.1021/ie401072x
% Problem 8

    % ---- Parámetros del problema ----
    cantidad_variables = 6;
    P = 40.52;
    T = 150.9;

    % Composición global (z)
    z = [0.5479, 0.0708, 0.0367, 0.0208, 0.0198, 0.304];
    %        C1      C2      C3     nC4     nC5      N2

    % Propiedades críticas y factor acéntrico (conjunto "Rosa", como en tu Fortran)
    Tc = [190.6, 305.4, 369.8, 425.2, 469.7, 126.2];
    Pc = [46.00, 48.80, 42.50, 38.00, 33.70, 33.90];
    w  = [0.008, 0.098, 0.153, 0.199, 0.251, 0.039];

    % Matriz de interacción binaria k_ij
    k = [ ...
        0.00, 0.00, 0.00, 0.00, 0.00, 0.02; ...
        0.00, 0.00, 0.00, 0.00, 0.00, 0.06; ...
        0.00, 0.00, 0.00, 0.00, 0.00, 0.08; ...
        0.00, 0.00, 0.00, 0.00, 0.00, 0.08; ...
        0.00, 0.00, 0.00, 0.00, 0.00, 0.08; ...
        0.02, 0.06, 0.08, 0.08, 0.08, 0.00 ];

    % ---- Construcción de x a partir de y ----
    % 5 variables libres: y(1:5). La sexta se cierra por suma = 1.
    x = zeros(1, cantidad_variables);
    x(1:cantidad_variables-1) = [y(1), y(2), y(3), y(4), y(5)];
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
        % Nota: se mantiene fidelidad al Fortran; x == 0 no se penaliza explícitamente.
        respuesta = sum( x .* ( log(x) + ln_phi_x - log(z) - ln_phi_z ) );
    end
end
