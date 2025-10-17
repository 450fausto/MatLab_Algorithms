function ln_phi = ln_phi_srk(x, Tc, Pc, w, k, P, T)
% ln_phi_srk  Fugacidades (ln phi_i) con la EoS SRK.
% Traducción fiel del Fortran dado.
%
% Entradas:
%   x  : fracciones molares (1xn o nx1)
%   Tc : T crítica (1xn)
%   Pc : P crítica (1xn)
%   w  : factor acéntrico (1xn)
%   k  : matriz k_ij (nxn)
%   P  : presión (escalar)
%   T  : temperatura (escalar)
%
% Salida:
%   ln_phi : (1xn) con ln(phi_i)

    % --- Asegurar vectores fila para operaciones vectorizadas
    if size(x,1)  > 1, x  = x.';  end
    if size(Tc,1) > 1, Tc = Tc.'; end
    if size(Pc,1) > 1, Pc = Pc.'; end
    if size(w,1)  > 1, w  = w.';  end

    n = numel(x);

    % --- Parámetros alfa, A, B (SRK)
    % alpha_i = [1 + (0.480 + 1.574 w_i - 0.176 w_i^2)*(1 - sqrt(T/Tc_i))]^2
    alpha = ( 1.0 + (0.480 + 1.574.*w - 0.176.*w.^2) .* (1 - sqrt(T./Tc)) ).^2;

    % A_i = 0.42747 * alpha_i * ( P * Tc_i^2 / (Pc_i * T^2) )
    A = 0.42747 .* alpha .* ( P .* (Tc.^2) ./ (Pc .* (T.^2)) );

    % B_i = 0.08664 * ( P * Tc_i / (Pc_i * T) )
    B = 0.08664 .* ( P .* Tc ./ (Pc .* T) );

    % Mezclas
    Bmix = sum(x .* B);

    % Am(i,j) = (1 - k_ij) * sqrt(A_i A_j)
    % Amix = sum_i sum_j x_i x_j Am_ij = x * Am * x'
    sqrtA = sqrt(A(:));                 % columna
    Am = (1 - k) .* (sqrtA * sqrtA.');  % nxn
    Amix = x * Am * x.';                % escalar

    % --- Ecuación cúbica de compresibilidad SRK:
    % z^3 - z^2 + (Amix - Bmix - Bmix^2) z - Amix*Bmix = 0
    termino_a = 1;
    termino_b = -1;
    termino_c = Amix - Bmix - Bmix.^2;
    termino_d = -Amix .* Bmix;

    % Raíces (Cardano) y elección de z (aquí, mínimo real como en Fortran: minval)
    raices = cardano(termino_a, termino_b, termino_c, termino_d);
    z = min(raices);  % se mantiene la selección del mínimo

    % --- Partes de ln(phi)
    % parte_a_i = (B_i/Bmix)*(z - 1) - ln(z - Bmix)
    parte_a = (B ./ Bmix) * (z - 1) - log(z - Bmix);

    % parte_b_i = (Amix/Bmix)*((B_i/Bmix) - 2 * sum_j Am_{i,j} x_j / Amix) * ln(1 + Bmix/z)
    % vectorizamos la suma fila i de Am(i,:).*x
    S = (Am * x.').';                    % S(i) = sum_j Am(i,j)*x_j  (fila a escalar)
    factor = (Amix / Bmix) .* ( (B./Bmix) - 2 * (S ./ Amix) );
    parte_b = factor * log(1 + Bmix / z);

    ln_phi = parte_a + parte_b;

    % Devolver orientación consistente con x
    if size(x,1) > 1, ln_phi = ln_phi.'; end
end
