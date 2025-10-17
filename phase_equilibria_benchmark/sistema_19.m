function respuesta = sistema_19(y)
% Agua + CO2 + 2-propanol + etanol, SRK
% x* = 0.18504982605063439, 0.0023796701744409404, 0.45441226048389355
% f* = -0.012649603564122792
%
% Harding, S. T., & Floudas, C. A. (2000).
% Phase stability with cubic equations of state: Global optimization approach.
% AIChE Journal, 46(7), 1422–1440. doi:10.1002/aic.690460715
% Example 3
%
% https://www.kaylaiacovino.com/Petrology_Tools/Critical_Constants_and_Acentric_Factors.htm

    % ---- Parámetros del problema ----
    cantidad_variables = 4;
    P = 22.5;
    T = 350;

    % Composición global (z)
    z = [0.99758, 0.00003, 0.00013, 0.00226];

    % Constantes críticas y factor acéntrico
    Tc = [647.3, 304.1, 508.3, 513.9];
    Pc = [221.2,  73.8,  47.6,  61.4 ];
    w  = [0.344,  0.239,  0.665, 0.644];

    % Interacciones binarias (aquí todas 0)
    k = zeros(cantidad_variables);

    % ---- Construcción de x a partir de y ----
    % 3 variables libres: y(1:3). La cuarta se cierra por suma = 1.
    x = [ y(1), y(2), y(3), 1 - y(1) - y(2) - y(3) ];

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
