function respuesta = sistema_36(y)
% N2 + C2, PR 
% x* = 0.48932300076583929
% f* = -0.013770597285383166
%
% Hua, J. Z., Brennecke, J. F., & Stadtherr, M. A. (1998).
% Reliable computation of phase stability using interval analysis.
% Computers & Chemical Engineering, 22(9), 1207–1214. doi:10.1016/s0098-1354(98)00024-6
% Problem 3 
%
% Ivanov, B. B., Galushko, A. A., & Stateva, R. P. (2013).
% Phase Stability Analysis with Equations of State—A Fresh Look from  a Different Perspective.
% Industrial & Engineering Chemistry Research, 52(32), 11208–11223. doi:10.1021/ie401072x
% Problem 3

    % ---- Parámetros del problema ----
    cantidad_variables = 2; %#ok<NASGU>
    P = 76.0; 
    T = 270.0;

    % Composición global (z)
    z = [0.3, 0.7];

    % Propiedades críticas y factor acéntrico
    Tc = [126.2, 305.4];
    Pc = [33.9,  48.8 ];
    w  = [0.04,  0.098];

    % Interacciones binarias k_ij
    k = [ ...
        0.0, 0.08; ...
        0.08, 0.0 ];

    % ---- Construcción de x a partir de y ----
    % Dos componentes: x = [x1, x2] con x2 = 1 - x1.
    x = [ y(1), 1 - y(1) ];

    % ---- ln(phi) con Peng–Robinson ----
    ln_phi_x = ln_phi_pr(x, Tc, Pc, w, k, P, T);
    ln_phi_z = ln_phi_pr(z, Tc, Pc, w, k, P, T);

    % Alinear dimensiones por si ln_phi_pr devuelve columna
    if size(ln_phi_x,1) > 1, ln_phi_x = ln_phi_x.'; end
    if size(ln_phi_z,1) > 1, ln_phi_z = ln_phi_z.'; end

    % ---- Penalización por fracciones negativas (emula PACK de Fortran) ----
    negativos = x(x < 0);
    criterio  = sum(negativos);   % < 0 si existe algún negativo

    % ---- Función objetivo de estabilidad (tangente de Gibbs) ----
    if criterio < 0
        respuesta = -criterio * 1.0e6;   % 10.0E5 en Fortran
    else
        % Nota: fidelidad al Fortran; x == 0 no se penaliza explícitamente.
        respuesta = sum( x .* ( log(x) + ln_phi_x - log(z) - ln_phi_z ) );
    end
end
