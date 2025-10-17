function respuesta = sistema_06(y)
% C1 + C2 + C3 + C4 + C5 + C6 + C7 + C8 + C9 + C10, SRK
% x* = 0.61410058753526386, 0.074469074863989237, 0.048200677764320229, 0.043016539266069562, 0.031783299830690830, 0.015123336264596665, 0.033972278680988968, 0.042539383539403126, 0.048644720385782388
% f* = -0.0000096609553336807350
%
% Bonilla-Petriciolet, A., Vázquez-Román, R., Iglesias-Silva, G. A., & Hall, K. R. (2006).
% Performance of Stochastic Global Optimization Methods in the Calculation of Phase Stability Analyses for Nonreactive and Reactive Mixtures.
% Industrial & Engineering Chemistry Research, 45(13), 4764–4772. doi:10.1021/ie051081g
%
% Ivanov, B. B., Galushko, A. A., & Stateva, R. P. (2013).
% Phase Stability Analysis with Equations of State—A Fresh Look from  a Different Perspective.
% Industrial & Engineering Chemistry Research, 52(32), 11208–11223. doi:10.1021/ie401072x
% Problem 11

    % ---- Parámetros del problema ----
    cantidad_variables = 10;
    P = 191.50;
    T = 435.35;

    % Composición global (z)
    z = [0.6436, 0.0752, 0.0474, 0.0412, 0.0297, 0.0138, 0.0303, 0.0371, 0.0415, 0.0402];
    %             C1      C2      C3      C4      C5      C6      C7      C8      C9      C10

    % Propiedades críticas y factor acéntrico (versión “precisa” del Fortran)
    Tc = [190.564, 305.322, 369.89, 425.125, 469.7000, 507.8200, 540.1300, 568.7400, 594.5500, 617.7000];
    Pc = [45.9920, 48.7220, 42.512, 37.9600, 33.67519, 30.44115, 27.36000, 24.83591, 22.81000, 21.03000];
    w  = [0.01142, 0.09900, 0.1521, 0.20081, 0.251032, 0.300319, 0.349000, 0.397528, 0.443300, 0.488400];

    % Rosa:
    % Tc = [190.6, 305.4, 369.8, 425.2, 469.7, 507.5, 540.3, 568.8, 594.6, 617.7];
    % Pc = [46.0,  48.80, 42.50, 38.00, 33.70, 30.10, 27.40, 24.90, 22.90, 21.20];
    % w  = [0.008, 0.098, 0.153, 0.199, 0.251, 0.299, 0.349, 0.398, 0.445, 0.489];

    % Interacción binaria k_ij (todo cero, como en tu Fortran)
    k = zeros(cantidad_variables);

    % ---- Construcción de x a partir de y ----
    % 9 variables libres: y(1:9). La décima se cierra por suma = 1.
    x = zeros(1, cantidad_variables);
    x(1:cantidad_variables-1) = [y(1), y(2), y(3), y(4), y(5), y(6), y(7), y(8), y(9)];
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
        respuesta = -criterio * 1.0e6;   % 10.0E5
    else
        % Nota: se mantiene fidelidad al Fortran; no se penaliza x == 0.
        respuesta = sum( x .* ( log(x) + ln_phi_x - log(z) - ln_phi_z ) );
    end
end
