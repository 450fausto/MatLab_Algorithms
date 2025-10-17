function respuesta = sistema_02(y)
% N2 + C1 + C2, SRK
% x* = 0.13682228197692597, 6.8392258985721613E-002
% f* = -0.014901225016232623
%
% Bonilla-Petriciolet, A., Vázquez-Román, R., Iglesias-Silva, G. A., & Hall, K. R. (2006).
% Performance of Stochastic Global Optimization Methods in the Calculation of Phase Stability Analyses for Nonreactive and Reactive Mixtures.
% Industrial & Engineering Chemistry Research, 45(13), 4764–4772. doi:10.1021/ie051081g
%
% Ivanov, B. B., Galushko, A. A., & Stateva, R. P. (2013).
% Phase Stability Analysis with Equations of State—A Fresh Look from a Different Perspective.
% Industrial & Engineering Chemistry Research, 52(32), 11208–11223. doi:10.1021/ie401072x
% Problem 5

    % ---- Parámetros del problema ----
    % Número de variables de decisión independientes (las demás se cierran con la restricción de suma)
    cantidad_variables = 3;

    % Condiciones de operación
    P = 76;      % Presión
    T = 270;     % Temperatura

    % Composición global (z), propiedades críticas (Tc, Pc), factor acéntrico (w)
    % y parámetros de interacción binaria k_ij (SRK)
    z  = [0.3,    0.1,     0.6   ];
    Tc = [126.2,  190.564, 305.32];
    Pc = [33.9,   45.9,    48.5  ];
    w  = [0.037,  0.011,   0.098 ];

    % Matriz simétrica de k_ij
    k = [ 0.000, 0.038, 0.080; ...
          0.038, 0.000, 0.021; ...
          0.080, 0.021, 0.000 ];

    % ---- Construcción de la composición tentativa x a partir de y ----
    % Para 3 componentes: x = [x1, x2, x3] con x3 = 1 - x1 - x2.
    x = [ y(1), y(2), 1 - y(1) - y(2) ];

    % ---- Cálculo de ln(phi) para x y para z usando SRK ----
    ln_phi_x = ln_phi_srk(x, Tc, Pc, w, k, P, T);
    ln_phi_z = ln_phi_srk(z, Tc, Pc, w, k, P, T);

    % Alinear dimensiones en fila para operaciones vectoriales
    if size(ln_phi_x,1) > 1, ln_phi_x = ln_phi_x.'; end
    if size(ln_phi_z,1) > 1, ln_phi_z = ln_phi_z.'; end

    % ---- Criterio de penalización por fracciones negativas ----
    % Emulamos PACK de Fortran: sumamos sólo los elementos negativos de x.
    negativos = x(x < 0);
    criterio  = sum(negativos);  % <= 0; es <0 si existe algún negativo

    % ---- Función objetivo de estabilidad (tangente de Gibbs) ----
    % Misma penalización que en Fortran: -criterio * 10.0E5 (= 1.0e6)
    if criterio < 0
        respuesta = -criterio * 1.0e6;
    else
        % Nota: se mantiene la fidelidad al Fortran (no se penaliza x == 0).
        % Si deseas evitar log(0), podemos usar max(x, eps) o penalizar x <= 0.
        respuesta = sum( x .* ( log(x) + ln_phi_x - log(z) - ln_phi_z ) );
    end

end
