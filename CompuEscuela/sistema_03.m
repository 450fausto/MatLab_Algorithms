function respuesta = sistema_03(y)
% El resultado debería ser -0.002 aprox
% C1 + C2 + C3 + C4 + C5 + C6 + C7–16 + C17+, SRK_2
% x* = 0.56458724554323736, 0.084643516682971098, 0.048310633773019709, 0.038244163849702834, 0.016974347676854337, 0.020305558881579006, 0.14615942932986148
% f* = -0.0062802362621926207
% SRK
% x* = 0.56789008804615304        8.4868666486083072E-002   4.8303528222366610E-002   3.8163889901578471E-002   1.6919196510887108E-002   2.0209367697312051E-002  0.14487354116496112
% f* = -5.7640611094150612E-003
%
% Bonilla-Petriciolet, A., Vázquez-Román, R., Iglesias-Silva, G. A., & Hall, K. R. (2006).
% Performance of Stochastic Global Optimization Methods in the Calculation of Phase Stability Analyses for Nonreactive and Reactive Mixtures.
% Industrial & Engineering Chemistry Research, 45(13), 4764–4772. doi:10.1021/ie051081g
%
% Ivanov, B. B., Galushko, A. A., & Stateva, R. P. (2013).
% Phase Stability Analysis with Equations of State—A Fresh Look from  a Different Perspective.
% Industrial & Engineering Chemistry Research, 52(32), 11208–11223. doi:10.1021/ie401072x
% Problem 9
%
% Harding, S. T., & Floudas, C. A. (2000).
% Phase stability with cubic equations of state: Global optimization approach.
% AIChE Journal, 46(7), 1422–1440. doi:10.1002/aic.690460715
% Example 4

    % ---- Parámetros generales ----
    cantidad_variables = 8;

    % Condiciones de operación
    P = 385.0;     % Presión
    T = 353.0;     % Temperatura

    % Composición global (z)
    z  = [0.7212, 0.09205, 0.04455, 0.03123, 0.01273, 0.01361, 0.07215, 0.01248];

    % Propiedades críticas y factor acéntrico
    Tc = [190.6, 305.4, 369.8, 425.2, 469.7, 507.5, 606.28, 835.67];
    Pc = [46.00, 48.80, 42.50, 38.00, 33.70, 30.10,  26.10,  14.77];
    w  = [0.008, 0.098, 0.153, 0.199, 0.251, 0.299,  0.4019, 0.7987];

    % Matriz de interacciones binarias k_ij (SRK)
    k = [ ...
        0.000, 0.002, 0.017, 0.015, 0.020, 0.039, 0.050, 0.090; ...
        0.002, 0.000, 0.000, 0.025, 0.010, 0.056, 0.040, 0.055; ...
        0.017, 0.000, 0.000, 0.000, 0.000, 0.000, 0.010, 0.010; ...
        0.015, 0.025, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000; ...
        0.020, 0.010, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000; ...
        0.039, 0.056, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000; ...
        0.050, 0.040, 0.010, 0.000, 0.000, 0.000, 0.000, 0.000; ...
        0.090, 0.055, 0.010, 0.000, 0.000, 0.000, 0.000, 0.000  ];

    % ---- Construcción de x a partir de y ----
    % 7 variables libres: y(1:7). La octava se cierra por suma = 1.
    x = zeros(1, cantidad_variables);
    x(1:cantidad_variables-1) = [y(1), y(2), y(3), y(4), y(5), y(6), y(7)];
    x(cantidad_variables)     = 1 - sum(x);

    % ---- Cálculo de ln(phi) con SRK ----
    ln_phi_x = ln_phi_srk(x, Tc, Pc, w, k, P, T);
    ln_phi_z = ln_phi_srk(z, Tc, Pc, w, k, P, T);

    % Alinear dimensiones (fila) para operaciones vectoriales
    if size(ln_phi_x,1) > 1, ln_phi_x = ln_phi_x.'; end
    if size(ln_phi_z,1) > 1, ln_phi_z = ln_phi_z.'; end

    % ---- Criterio de penalización por fracciones negativas ----
    % Emulación de PACK de Fortran: sumar sólo los elementos negativos de x.
    negativos = x(x < 0);
    criterio  = sum(negativos);   % <= 0; es <0 si existe algún negativo

    % ---- Función objetivo de estabilidad (tangente de Gibbs) ----
    % Misma penalización que en Fortran: -criterio * 10.0E5 (= 1.0e6)
    if criterio < 0
        respuesta = -criterio * 1.0e6;
    else
        % Nota: se mantiene fidelidad al Fortran; no se penaliza x == 0.
        % Si quieres evitar log(0), podemos usar max(x, realmin) o penalizar x <= 0.
        respuesta = sum( x .* ( log(x) + ln_phi_x - log(z) - ln_phi_z ) );
    end

end
