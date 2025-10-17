function respuesta = sistema_20(w)
% Ácido acético + benceno + furfural + ciclohexano, UNIQUAC
% x* =  0.017514477936952964, 0.19953729228166880, 0.13392662040944914
% f* = -0.0049317723847621814
%
% Tessier S. R., Brennecke J. F. and Stadtherr M. A. (2000).
% Reliable phase stability analysis for excess Gibbs energy models.
% Chemical Engineering Science, 55(10), 1785-1796. Doi: 10.1016/s0009-2509(99)00442-x
% Problem 5

    % ---- Parámetros del problema ----
    cantidad_variables = 4;

    % Composición global (z)
    z = [0.05, 0.20, 0.35, 0.40];

    % Parámetros UNIQUAC
    r  = [2.2024, 3.1878, 3.1680, 4.0464];
    q  = [2.072,  2.400,  2.484,  3.240 ];
    qq = [2.072,  2.400,  2.484,  3.240 ];

    % Matriz tau (tal como en el Fortran)
    tau = [ ...
        1.00000, 1.26362, 3.36860, 0.85128; ...
        0.99972, 1.00000, 1.02041, 0.89333; ...
        0.31633, 0.79027, 1.00000, 0.96249; ...
        0.49739, 1.09619, 0.26222, 1.00000 ];

    % ---- Construcción de x a partir de w ----
    % Cuatro componentes: x = [x1, x2, x3, x4] con x4 = 1 - x1 - x2 - x3.
    x = [ w(1), w(2), w(3), 1 - w(1) - w(2) - w(3) ];

    % ---- ln(gamma) con UNIQUAC ----
    lngx = ln_gamma_uniquac(x, r, q, qq, tau);
    lngz = ln_gamma_uniquac(z, r, q, qq, tau);

    % Alinear dimensiones a fila por si ln_gamma_uniquac devuelve columna
    if size(lngx,1) > 1, lngx = lngx.'; end
    if size(lngz,1) > 1, lngz = lngz.'; end

    % ---- Penalización por fracciones negativas ----
    negativos = x(x < 0);
    criterio  = sum(negativos);   % <= 0; <0 si existe algún negativo

    % ---- Función objetivo de estabilidad (tangente de Gibbs) ----
    if criterio < 0
        respuesta = -criterio * 1.0e6;   % 10.0E5 en Fortran
    else
        % Nota: fidelidad al Fortran; x == 0 no se penaliza explícitamente.
        respuesta = sum( x .* ( log(x) + lngx - log(z) - lngz ) );
    end
end
