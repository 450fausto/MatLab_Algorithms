function ln_phi = ln_phi_srk_2(x, Tc, Pc, w, k, P, T)
% ln_phi_srk_2  Fugacidades (ln phi_i) con SRK usando las correlaciones
% de componente puro de Graboski–Daubert para alpha(T, w).
%
% Entradas:
%   x  : fracciones molares (1xn o nx1)
%   Tc : temperatura crítica (1xn)
%   Pc : presión crítica (1xn)
%   w  : factor acéntrico (1xn)
%   k  : matriz k_ij (nxn)
%   P  : presión (escalar)
%   T  : temperatura (escalar)
%
% Salida:
%   ln_phi : (1xn) con ln(phi_i)
%
% Nota: aquí se selecciona la **máxima** raíz real de z (fase vapor),
% igual que en tu Fortran (maxval).

    % --- Alinear orientación a fila
    if size(x,1)  > 1, x  = x.';  end
    if size(Tc,1) > 1, Tc = Tc.'; end
    if size(Pc,1) > 1, Pc = Pc.'; end
    if size(w,1)  > 1, w  = w.';  end

    n = numel(x);

    % --- Parámetros SRK (Graboski–Daubert para alpha)
    % alpha_i = [1 + (0.48508 + 1.55171 w_i - 0.15613 w_i^2)*(1 - sqrt(T/Tc_i))]^2
    alpha = ( 1.0 + (0.48508 + 1.55171.*w - 0.15613.*w.^2) .* (1 - sqrt(T./Tc)) ).^2;

    % A_i = 0.42747 * alpha_i * ( P * Tc_i^2 / (Pc_i * T^2) )
    A = 0.42747 .* alpha .* ( P .* (Tc.^2) ./ (Pc .* (T.^2)) );

    % B_i = 0.08664 * ( P * Tc_i / (Pc_i * T) )
    B = 0.08664 .* ( P .* Tc ./ (Pc .* T) );

    % --- Mezclas
    Bmix = sum(x .* B);

    % Am(i,j) = (1 - k_ij) * sqrt(A_i A_j)
    sqrtA = sqrt(A(:));                  % columna
    Am = (1 - k) .* (sqrtA * sqrtA.');   % nxn
    Amix = x * Am * x.';                 % escalar

    % --- Cúbica de SRK para z:
    % z^3 - z^2 + (Amix - Bmix - Bmix^2) z - Amix*Bmix = 0
    a1 = 1;
    a2 = -1;
    a3 = Amix - Bmix - Bmix.^2;
    a4 = -Amix .* Bmix;

    raices = cardano(a1, a2, a3, a4);
    z = max(raices);   % Fase vapor (como en Fortran: maxval)

    % --- ln(phi)
    % parte_a_i = (B_i/Bmix)*(z - 1) - ln(z - Bmix)
    parte_a = (B ./ Bmix) * (z - 1) - log(z - Bmix);

    % parte_b_i = (Amix/Bmix)*((B_i/Bmix) - 2 * sum_j Am_{i,j} x_j / Amix) * ln(1 + Bmix/z)
    S = (Am * x.').';                         % S(i) = sum_j Am(i,j)*x_j
    factor = (Amix / Bmix) .* ( (B./Bmix) - 2*(S./Amix) );
    parte_b = factor * log(1 + Bmix / z);

    ln_phi = parte_a + parte_b;

    % Regresar como columna si x entró columna
    if size(x,1) > 1, ln_phi = ln_phi.'; end
end
