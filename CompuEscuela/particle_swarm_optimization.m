function [solucion,convergencia]=particle_swarm_optimization(funcion_costo,limites_inferiores,limites_superiores,tamanio_poblacion,maximo_iteraciones)
%% Parámetros de control
w = 1;            % Inertia Weight
wdamp = 0.99;     % Inertia Weight Damping Ratio
c1 = 1.5;         % Personal Learning Coefficient
c2 = 2.0;         % Global Learning Coefficient

%% Matrices Iniciales
numero_variables = length(limites_superiores);
convergencia=zeros(maximo_iteraciones,1);
costo = zeros(1, tamanio_poblacion);
velocidad = zeros(tamanio_poblacion, numero_variables);
mejor_posicion = zeros(tamanio_poblacion, numero_variables);
mejor_costo = zeros(1, tamanio_poblacion);
posicion = repmat(limites_inferiores, tamanio_poblacion, 1) + ...
    (repmat(limites_superiores, tamanio_poblacion, 1) - ...
    repmat(limites_inferiores,tamanio_poblacion, 1))...
    .* rand(tamanio_poblacion, numero_variables);

%% Inicialización
for i=1:tamanio_poblacion
    costo(i)=funcion_costo(posicion(i,:));
    mejor_posicion(i,:) = posicion(i,:);
    mejor_costo(i) = costo(i);
end
[mejor_costo_global,b] = min(costo);
mejor_posicion_global = posicion(b,:);
convergencia(1) = mejor_costo_global;
velocidad_maxima = 0.1 * (limites_superiores - limites_inferiores);
velocidad_minima = -velocidad_maxima;

%% Proceso iterativo
for iter=1:maximo_iteraciones
    for pop=1:tamanio_poblacion
        velocidad(pop,:) = w*velocidad(pop,:) ...
            +c1*rand(1, numero_variables).*(mejor_posicion(pop) - posicion(pop,:)) ...
            +c2*rand(1, numero_variables).*(mejor_posicion_global - posicion(pop,:));
        velocidad(pop,:) = max(min(velocidad(pop,:), velocidad_maxima), velocidad_minima);
        posicion(pop,:) = posicion(pop,:) + velocidad(pop,:);
        esta_fuera = posicion(pop,:) > limites_superiores | posicion(pop,:) < limites_inferiores;
        velocidad(pop,esta_fuera) = -velocidad(pop,esta_fuera);
        posicion(pop,:) = max(min(posicion(pop,:), limites_superiores), limites_inferiores);
        costo(pop) = funcion_costo(posicion(pop));
        if costo(pop) < mejor_costo(i)
            mejor_posicion(pop,:) = posicion(pop,:);
            mejor_costo(pop) = costo(pop);
            if mejor_costo(pop) < mejor_costo_global
                mejor_posicion_global = mejor_posicion(pop,:);
                mejor_costo_global = mejor_costo(pop);
            end
        end
    end
    convergencia(iter) = mejor_costo_global;
    w = w*wdamp;
end
%% salida
solucion = [mejor_posicion_global, mejor_costo_global];
end