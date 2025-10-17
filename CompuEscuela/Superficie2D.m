% Curva 2D del funcional de estabilidad: sistema_12
% Dominio: y = [x1], con x2 = 1 - x1  (0 < x1 < 1)

clear; clc;

% --- Malla 1D en x1 ---
N = 2001;               % puntos (ajusta si quieres más/menos)
epsx = 1e-12;           % evita evaluar exactamente en 0 o 1
x1 = linspace(epsx, 1-epsx, N);

% --- Evaluar sistema_12 en cada x1 ---
F = arrayfun(@(u) sistema_27(u), x1);

% --- Graficar ---
figure('Name','Curva sistema_12','Color','w');
plot(x1, F, 'LineWidth', 1.5);
xlabel('x_1'); ylabel('f(x_1)');
title('Sistema 12: C1 + C3 (SRK)');
grid on; box on;

% --- Marcar el punto reportado (opcional) ---
x1_star = 0.0;
f_star  = sistema_27(x1_star);
hold on; plot(x1_star, f_star, 'r.', 'MarkerSize', 16);
legend({'f(x_1)','punto reportado'}, 'Location','best');
hold off;

% --- Encontrar y marcar el minimo numerico (opcional) ---
try
    fun = @(u) sistema_27(u);
    opts = optimset('TolX',1e-10,'Display','off');
    [xmin, fmin] = fminbnd(fun, epsx, 1-epsx, opts);
    hold on; plot(xmin, fmin, 'ko', 'MarkerSize', 6, 'MarkerFaceColor','y');
    text(xmin, fmin, sprintf('  min: x_1=%.6f, f=%.6g', xmin, fmin), 'VerticalAlignment','bottom');
    hold off;
catch
    % Si no tienes Optimization Toolbox, omite este bloque
end
