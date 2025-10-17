function respuesta = sistema_01(y)
% H2S + C1, SRK
% x* = 7.6619630769793529E-002
% f* = -3.9319962731680751E-003
%
% Hua, J. Z., Brennecke, J. F., & Stadtherr, M. A. (1998).
% Reliable computation of phase stability using interval analysis.
% Computers & Chemical Engineering, 22(9), 1207–1214. doi:10.1016/S0098-1354(98)00024-6
% Problem 1
%
% Ivanov, B. B., Galushko, A. A., & Stateva, R. P. (2013).
% Phase Stability Analysis with Equations of State—A Fresh Look from a Different Perspective.
% Industrial & Engineering Chemistry Research, 52(32), 11208–11223. doi:10.1021/ie401072x
% Problem 1
%
% Harding, S. T., & Floudas, C. A. (2000).
% Phase stability with cubic equations of state: Global optimization approach.
% AIChE Journal, 46(7), 1422–1440. doi:10.1002/aic.690460715
% Example 2

    % ---- Definición de variables del problema ----
    % Cantidad de variables de decisión (una fracción independiente; la otra es 1 - y(1))
    cantidad_variables = 2;

    % Condiciones de operación
    P = 40.53;   % Presión
    T = 190;     % Temperatura

    % Composición global (z), propiedades críticas (Tc, Pc), factor acéntrico (w)
    % y parámetros de interacción binaria (kij) para SRK
    z  = [0.0187, 0.9813];
    Tc = [373.2, 190.6];
    Pc = [89.4,   46.0 ];
    w  = [0.1,    0.008];

    % Matriz de k_ij (simétrica en este caso)
    k = [ 0.00, 0.08; ...
          0.08, 0.00 ];

    % ---- Construcción de la composición tentativa x a partir de y ----
    % En este caso, hay 2 componentes: x = [x1, x2] con x2 = 1 - x1.
    % y(1) es la variable libre.
    x = [y(1), 1 - y(1)];

    % ---- Cálculo de ln(phi) para x y para z usando SRK ----
    % ln_phi_srk debe devolver un vector 1x2 (o 2x1); aquí trabajamos como fila.
    ln_phi_x = ln_phi_srk(x, Tc, Pc, w, k, P, T);
    ln_phi_z = ln_phi_srk(z, Tc, Pc, w, k, P, T);

    % Asegurar que las dimensiones sean consistentes (filas) para operaciones elemento a elemento
    if size(ln_phi_x,1) > 1, ln_phi_x = ln_phi_x.'; end
    if size(ln_phi_z,1) > 1, ln_phi_z = ln_phi_z.'; end

    % ---- Criterio de penalización por fracciones negativas ----
    % Imitamos el PACK de Fortran: sumamos sólo los elementos negativos de x.
    negativos = x(x < 0);     % vector con las x_i negativas
    criterio  = sum(negativos); % será <= 0; <0 si hay algún negativo

    % ---- Función objetivo de estabilidad (Tangente de Gibbs) ----
    % Si hay negativos, penalizamos con 1.0E6 * (-criterio), igual que en el Fortran (10.0E5)
    if criterio < 0
        respuesta = -criterio * 1.0e6;
    else
        % Nota: en el código original no se penaliza x == 0, lo cual haría log(0) = -Inf.
        % Se mantiene el comportamiento de Fortran para fidelidad.
        respuesta = sum( x .* ( log(x) + ln_phi_x - log(z) - ln_phi_z ) );
    end

end
