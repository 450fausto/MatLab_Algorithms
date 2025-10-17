% Superficie 3D del funcional de estabilidad: sistema_07
% Dominio: w = [w1, w2], con x3 = 1 - w1 - w2 >= 0

clear; clc;

% --- Resolucion y dominio ---
N = 201;           % densidad de malla
epsx = 1e-12;      % evita evaluar exactamente en los bordes (log(0))

w1v = linspace(0, 1, N);
w2v = linspace(0, 1, N);
[W1, W2] = meshgrid(w1v, w2v);

% --- Region factible del simplex ---
W3 = 1 - W1 - W2;
mask = (W1 >= epsx) & (W2 >= epsx) & (W3 >= epsx);

% --- Evaluacion del funcional ---
F = NaN(size(W1));
F(mask) = arrayfun(@(a,b) sistema_39([a, b]), W1(mask), W2(mask));

% --- Valor minimo finito para trazar el borde del simplex ---
if any(~isnan(F(:)))
    zline = min(F(~isnan(F)));
else
    zline = 0;
end

% --- Superficie 3D ---
figure('Name','Superficie sistema_07 (NRTL)','Color','w');
surf(W1, W2, F, 'EdgeColor','none');
xlabel('x_1 (w_1)','Interpreter','none');
ylabel('x_2 (w_2)','Interpreter','none');
zlabel('f(w)','Interpreter','none');
title('Superficie 3D de sistema_07 (Tolueno + Agua + Anilina, NRTL)', 'Interpreter','none');
colorbar; view(45, 30); grid on; box on; hold on;
plot3([0 1 0 0], [0 0 1 0], [zline zline zline zline], 'k-', 'LineWidth', 1);
hold off;

% --- Contornos 2D ---
figure('Name','Contornos sistema_07','Color','w');
contourf(W1, W2, F, 30, 'LineStyle','none'); colorbar;
xlabel('x_1 (w_1)','Interpreter','none');
ylabel('x_2 (w_2)','Interpreter','none');
title('Contornos de f(w) en el simplex', 'Interpreter','none');
axis equal tight; hold on;
plot([0 1 0 0], [0 0 1 0], 'k-', 'LineWidth', 1);
hold off;

% --- Punto reportado en el articulo/Fortran (opcional) ---
x1_star = 0.000066937752536364210;
x2_star = 0.99686525905361334;
f_star  = sistema_39([x1_star, x2_star]);

figure(findobj('Name','Contornos sistema_07'));
hold on; plot(x1_star, x2_star, 'r.', 'MarkerSize', 18); hold off;

figure(findobj('Name','Superficie sistema_07 (NRTL)'));
hold on; plot3(x1_star, x2_star, f_star, 'r.', 'MarkerSize', 18); hold off;
