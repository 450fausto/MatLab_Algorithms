function respuesta = sistema_12(y)
% C1 + C3, SRK
% x* = 0.77244380013981195
% f* = -0.00033438415311950467
%
% Hua, J. Z., Brennecke, J. F., & Stadtherr, M. A. (1998).
% Reliable computation of phase stability using interval analysis.
% Computers & Chemical Engineering, 22(9), 1207–1214. doi:10.1016/s0098-1354(98)00024-6
% Problem 2
%
% Ivanov, B. B., Galushko, A. A., & Stateva, R. P. (2013).
% Phase Stability Analysis with Equations of State—A Fresh Look from  a Different Perspective.
% Industrial & Engineering Chemistry Research, 52(32), 11208–11223. doi:10.1021/ie401072x
% Problem 2

    % ---- Parámetros del problema ----
    cantidad_variables = 2;
    P = 100.0;
    T = 277.6;

    % Composición global (z)
    z = [0.68, 0.32];

    % Propiedades críticas y factor acéntrico
    Tc = [190.6, 369.8];
    Pc = [46.0,  42.5 ];
    w  = [0.008, 0.152];

    % Matriz de interacción binaria k_ij
    k = [ 0.000, 0.029; ...
          0.029, 0.000 ];

    % ---- Construcción de x a partir de y ----
    % Dos componentes: x = [x1, x2] con x2 = 1 - x1.
    x = [ y(1), 1 - y(1) ];

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
        respuesta = -criterio * 1.0e6;   % 10.0E5 en Fortran
    else
        % Nota: se mantiene fidelidad al Fortran; x == 0 no se penaliza explícitamente.
        respuesta = sum( x .* ( log(x) + ln_phi_x - log(z) - ln_phi_z ) );
    end
end
