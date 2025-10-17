function respuesta = sistema_17(w)
% Ácido acético + benceno + furfural + ciclohexano + agua, UNIQUAC
% x* =  0.22748543782080269, 0.0033447179085991561, 0.047502524800731014, 0.0015842688961442879
% f* = -0.17765416298801889
%
% Tessier S. R., Brennecke J. F. and Stadtherr M. A. (2000).
% Reliable phase stability analysis for excess Gibbs energy models.
% Chemical Engineering Science, 55(10), 1785-1796. Doi: 10.1016/s0009-2509(99)00442-x
% Problem 6

    % ---- Parámetros del problema ----
    cantidad_variables = 5;

    % Composición global (z)
    z = [0.2, 0.2, 0.2, 0.2, 0.2];

    % Parámetros UNIQUAC
    r  = [2.2024, 3.1878, 3.1680, 4.0464, 0.9200];
    q  = [2.072,  2.400,  2.484,  3.240,  1.400 ];
    qq = [2.072,  2.400,  2.484,  3.240,  1.400 ];

    % Matriz tau (tal como en el Fortran)
    tau = [ ...
        1.00000, 1.26362, 3.36860, 0.85128, 1.54662; ...
        0.99972, 1.00000, 1.02041, 0.89333, 0.09441; ...
        0.31633, 0.79027, 1.00000, 0.96249, 0.60488; ...
        0.49739, 1.09619, 0.26222, 1.00000, 0.08839; ...
        2.44225, 0.13507, 0.69066, 0.19491, 1.00000 ];

    % ---- Construcción de x a partir de w ----
    % Cinco componentes: x = [x1, x2, x3, x4] y x5 = 1 - x1 - x2 - x3 - x4.
    x = [ w(1), w(2), w(3), w(4), 1 - w(1) - w(2) - w(3) - w(4) ];

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
