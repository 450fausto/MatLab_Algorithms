function raices = cardano(a, b, c, d)
% CARDANO  Resuelve ax^3 + bx^2 + cx + d = 0 con el m�todo de Cardano.
% Correcci�n: en el caso |delta| ~ 0 se aplica el corrimiento -b/(3a)
% (ra�z triple y/o doble+simple). Se a�ade clip num�rico al argumento de acos.

    pi_val = 3.141592653589793;
    raices = zeros(1, 3);

    % Par�metros de Cardano (forma deprimida t^3 + p t + q)
    p = (3*a*c - b^2) / (3*a^2);
    q = (2*b^3 - 9*a*b*c + 27*a^2*d) / (27*a^3);

    % Discriminante escalado (4?): signo equivalente al ? cl�sico
    delta = q^2 + (4 * p^3) / 27;

    shift = b / (3*a);  % corrimiento a x = t - b/(3a)

    if delta > 0
        % �nica ra�z real
        z1 = (-q + sqrt(delta)) / 2;
        z2 = (-q - sqrt(delta)) / 2;
        u  = crts(z1);
        v  = crts(z2);
        root = u + v - shift;
        raices(:) = root;  % mismo comportamiento que el Fortran original

    elseif abs(delta) < 1.0e-12
        % Caso degenerado (ra�z m�ltiple)
        if abs(p) < 1.0e-15
            % p = 0, q = 0  -> ra�z triple
            rt = -shift;
            raices(:) = rt;
        else
            % Una ra�z simple y una doble (en variable t), luego aplicar corrimiento
            t1 =  3*q / p;
            t2 = -3*q / (2*p);
            raices(1) = t1 - shift;
            raices(2) = t2 - shift;
            raices(3) = t2 - shift;
        end

    else
        % Tres ra�ces reales (f�rmula trigonom�trica)
        % Clip por estabilidad num�rica:
        arg = (3*q/(2*p)) * sqrt(-3/p);
        arg = max(-1, min(1, arg));

        tetha = acos(arg);
        r = 2 * sqrt(-p/3);

        raices(1) =  r * cos(           tetha/3)            - shift;
        raices(2) =  r * cos((tetha + 2*pi_val)/3)          - shift;
        raices(3) =  r * cos((tetha + 4*pi_val)/3)          - shift;
    end
end

% ---------- Auxiliar: ra�z c�bica real con signo ----------
function y = crts(z)
    if z == 0
        y = 0;
    else
        y = sign(z) * abs(z)^(1/3);
    end
end
