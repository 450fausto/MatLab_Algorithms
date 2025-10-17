function ln_gamma = ln_gamma_uniquac(x, r, q, qq, tau)
% ln_gamma_uniquac  (UNIQUAC) 
% Traducción directa del Fortran proporcionado, conservando comentarios y estructura.
%
% Entradas:
%   x   : vector fracciones molares (1xn o nx1)
%   r,q : parámetros UNIQUAC
%   qq  : parámetro adicional (q') usado en los términos residual/estructural de esta variante
%   tau : matriz (n x n) de interacciones (tau_ij)
%
% Salida:
%   ln_gamma : vector (misma forma que x) con ln(gamma_i)
%
% NOTA: se asume que log() = log natural (igual que en Fortran).
%       Se preserva el uso de qq en 'tetha_p' y en los términos D/E.
%
% ADVERTENCIA (del código original):
% cuando se usa shape en tau declarada como parameter, se convierte en la matriz
% transpuesta de tau, por esta razón se usan los índices al revés. ej. en lugar de tau_ij,
% se tiene que usar tau_ji

    tau = tau';

    % Asegurar vectores fila para operaciones vectorizadas
    if size(x,1) > 1, x = x.'; end
    if size(r,1) > 1, r = r.'; end
    if size(q,1) > 1, q = q.'; end
    if size(qq,1) > 1, qq = qq.'; end

    % ---- Parámetros y memoria ----
    z = 10;  % coordinación efectiva en UNIQUAC
    n = numel(x);

    ln_gamma  = zeros(1, n);
    phi       = zeros(1, n);
    tetha     = zeros(1, n);
    tetha_p   = zeros(1, n);
    L         = zeros(1, n);
    termino_a = zeros(1, n);
    termino_b = zeros(1, n);
    termino_c = zeros(1, n);
    termino_d = zeros(1, n);
    termino_e = zeros(1, n);

    % ---- Definiciones volumétricas y superficiales ----
    % phi_i    = r_i x_i / sum_j (r_j x_j)
    % tetha_i  = q_i x_i / sum_j (q_j x_j)
    % tetha'_i = qq_i x_i / sum_j (qq_j x_j)
    phi     = (r .* x) ./ sum(x .* r);
    tetha   = (q .* x) ./ sum(q .* x);
    tetha_p = (qq .* x) ./ sum(qq .* x);

    % L_i = (z/2)*(r_i - q_i) - (r_i - 1)
    L = (z/2) .* (r - q) - (r - 1);

    % ---- Términos combinatorio y estructural (según tu Fortran) ----
    % termino_a = log(phi ./ x)
    % termino_b = (z/2) * q .* log(tetha ./ phi)
    % termino_c = -(phi ./ x) * sum(x .* L)
    termino_a = log(phi ./ x);
    termino_b = (z/2) .* q .* log(tetha ./ phi);
    termino_c = - (phi ./ x) * sum(x .* L);

    % ---- Términos residuales con qq y tau ----
    % Ojo con los índices "al revés" como en tu nota.
    % termino_d(i) = - qq(i) * log( sum_j tetha'_j * tau(i,j) )
    % termino_e(i) = - qq(i) * sum_j [ tetha'_j * tau(j,i) / sum_k tetha'_k * tau(j,k) ]
    for i = 1:n
        % sum_j tetha_p(j) * tau(i,j)  (producto con FILA i de tau)
        s_row_i = sum(tetha_p .* tau(i, :));
        termino_d(i) = - qq(i) * log(s_row_i);

        % suma_j tetha_p(j) * tau(j,i) / sum_k tetha_p(k) * tau(j,k)
        acc = 0.0;
        for j = 1:n
            denom_j = sum(tetha_p .* tau(j, :));
            acc = acc + tetha_p(j) * tau(j, i) / denom_j;
        end
        termino_e(i) = - qq(i) * acc;
    end

    % ---- Ensamble final (tal cual en el Fortran) ----
    % ln_gamma = termino_a + termino_b + L + termino_c + termino_d + qq + termino_e
    ln_gamma = termino_a + termino_b + L + termino_c + termino_d + qq + termino_e;

    % Restaurar orientación a columna si x vino como columna
    if size(x,1) > 1, ln_gamma = ln_gamma.'; end
end
