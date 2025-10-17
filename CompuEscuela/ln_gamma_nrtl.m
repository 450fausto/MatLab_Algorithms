function ln_gamma = ln_gamma_nrtl(x, tau, G)
% ln_gamma_nrtl
% Traducción directa del Fortran dado (NRTL).
%
% Entradas:
%   x   : fracciones molares (1xn o nx1)
%   tau : matriz (n x n) de parámetros ?_ij
%   G   : matriz (n x n) de G_ij
%
% Salida:
%   ln_gamma : vector (1xn) con ln(gamma_i)
%
% Nota: en el Fortran aparece un parámetro z=10 que no se usa en esta rutina.
    % Parche, debido a que en Fortran se toma la matriz transpuesta
    tau = tau';
    G = G';
    % Asegurar vectores fila para operaciones vectorizadas
    
    if size(x,1) > 1, x = x.'; end

    n = numel(x);
    ln_gamma = zeros(1, n);

    for i = 1:n
        suma = 0.0;
        for j = 1:n
            % sumb = sum_k G_{j,k} * x_k
            sumb = sum( G(j, :) .* x );

            % sumc = sum_k tau_{j,k} * G_{j,k} * x_k
            sumc = sum( tau(j, :) .* G(j, :) .* x );

            % suma += x_j * ( tau_{j,i}*G_{j,i}*sumb - G_{j,i}*sumc ) / sumb^2
            % (respetando exactamente la forma del Fortran)
            suma = suma + x(j) * ( tau(j, i) * G(j, i) * sumb - G(j, i) * sumc ) / (sumb.^2);
        end

        % ln_gamma(i) = [sum_k tau_{i,k} G_{i,k} x_k] / [sum_k G_{i,k} x_k] + suma
        ln_gamma(i) = sum( tau(i, :) .* G(i, :) .* x ) / sum( G(i, :) .* x ) + suma;
    end

    % Conservar orientación: si x llegó como columna, devolver columna
    if size(x,1) > 1, ln_gamma = ln_gamma.'; end
end
