function respuesta = sistema_09(y)
% C2 + C3 + nC4 + nC5 + nC6, SRK
% x* = 0.38829635778801025, 0.29247235943977556, 0.20476729387874437, 0.074901399711355093
% f* = -0.0000020328388097537212
%
% Bonilla-Petriciolet, A., Vázquez-Román, R., Iglesias-Silva, G. A., & Hall, K. R. (2006).
% Performance of Stochastic Global Optimization Methods in the Calculation of Phase Stability Analyses for Nonreactive and Reactive Mixtures.
% Industrial & Engineering Chemistry Research, 45(13), 4764–4772. doi:10.1021/ie051081g

    % ---- Parámetros del problema ----
    cantidad_variables = 5;
    P = 55.83;
    T = 390.0;

    % Composición global (z)
    z = [0.401, 0.293, 0.199, 0.0707, 0.0363];
    %          C2     C3     C4     C5      C6

    % Propiedades críticas y factor acéntrico (versión precisa del Fortran)
    Tc = [305.322, 369.89, 425.125, 469.7000, 507.8200];
    Pc = [48.7220, 42.512, 37.9600, 33.67519, 30.44115];
    w  = [0.09900, 0.1521, 0.20081, 0.251032, 0.300319];

    % Interacción binaria k_ij (todo cero, como en el Fortran)
    k = zeros(cantidad_variables);

    % ---- Construcción de x a partir de y ----
    % 4 variables libres: y(1:4). La quinta se cierra por suma = 1.
    x = zeros(1, cantidad_variables);
    x(1:cantidad_variables-1) = [y(1), y(2), y(3), y(4)];
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
