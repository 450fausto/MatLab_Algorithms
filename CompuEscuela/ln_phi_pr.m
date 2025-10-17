function ln_phi = ln_phi_pr(x, Tc, Pc, w, k, P, T)
% ln_phi_pr  Calcula ln(phi_i) con la EoS de Peng–Robinson (PR).
% Traducción fiel del Fortran proporcionado.
%
% Entradas:
%   x  : fracciones molares (1xn o nx1)
%   Tc : temperatura crítica (1xn)
%   Pc : presión crítica (1xn)
%   w  : factor acéntrico (1xn)
%   k  : matriz de interacción binaria k_ij (nxn)
%   P  : presión (escalar)
%   T  : temperatura (escalar)
%
% Salida:
%   ln_phi : (1xn) con ln(phi_i)
%
% Nota: se selecciona la **mínima** raíz real de la cúbica (fase líquida),
% replicando el uso de minval(raices) del Fortran.

    % --- Alinear a vectores fila
    if size(x,1)  > 1, x  = x.';  end
    if size(Tc,1) > 1, Tc = Tc.'; end
    if size(Pc,1) > 1, Pc = Pc.'; end
    if size(w,1)  > 1, w  = w.';  end

    n = numel(x);

    % --- Parámetros PR (kappa de Soave para PR)
    % alpha_i = [1 + (0.37464 + 1.54226 w_i - 0.26992 w_i^2)*(1 - sqrt(T/Tc_i))]^2
    alpha = ( 1.0 + (0.37464 + 1.54226.*w - 0.26992.*w.^2) .* (1 - sqrt(T./Tc)) ).^2;

    % A_i = 0.45724 * alpha_i * ( P * Tc_i^2 / (Pc_i * T^2) )
    A = 0.45724 .* alpha .* ( P .* (Tc.^2) ./ (Pc .* (T.^2)) );

    % B_i = 0.0778 * ( P * Tc_i / (Pc_i * T) )
    B = 0.0778 .* ( P .* Tc ./ (Pc .* T) );

    % --- Mezclas
    Bmix = sum(x .* B);

    % Am(i,j) = (1 - k_ij) * sqrt(A_i A_j)
    sqrtA = sqrt(A(:));                    % columna
    Am = (1 - k) .* (sqrtA * sqrtA.');     % nxn
    Amix = x * Am * x.';                   % escalar

    % --- Cúbica de PR (tal cual tu Fortran):
    % z^3 + (Bmix - 1) z^2 + (Amix - 2 Bmix - 3 Bmix^2) z + (Bmix^3 + Bmix^2 - Amix Bmix) = 0
    a1 = 1;
    a2 = (Bmix - 1);
    a3 = (Amix - 2*Bmix - 3*Bmix.^2);
    a4 = (Bmix.^3 + Bmix.^2 - Amix.*Bmix);

    raices = cardano(a1, a2, a3, a4);
    z = min(raices);   % fase líquida (minval)

    % --- ln(phi)
    % parte_a_i = (B_i/Bmix)*(z - 1) - ln(z - Bmix)
    parte_a = (B ./ Bmix) * (z - 1) - log(z - Bmix);

    % parte_b_i = (Amix/Bmix) * ((B_i/Bmix) - 2 * sum_j Am_{i,j} x_j / Amix)
    %             * ln( (z + (1+sqrt(2)) Bmix) / (z + (1-sqrt(2)) Bmix) )
    S = (Am * x.').';                                     % S(i) = sum_j Am(i,j)*x_j
    log_arg = ( z + (1 + sqrt(2)) * Bmix ) / ( z + (1 - sqrt(2)) * Bmix );
    factor = (Amix / Bmix) .* ( (B./Bmix) - 2*(S./Amix) );
    parte_b = factor * log(log_arg);

    ln_phi = parte_a + parte_b / sqrt(8.0);

    % Devolver como columna si x entró columna
    if size(x,1) > 1, ln_phi = ln_phi.'; end
end
