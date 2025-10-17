function [solucion, convergencia] = exam_copying_room(funcion_costo, limites_inferiores, limites_superiores, tamanio_poblacion, maximo_iteraciones)
%% ---------- Parámetros de control ----------
temperatura_social = 0.7;                % temperatura para softmax de consejos
tasa_confianza = 0.15;                 % tasa de aprendizaje de confianza
tolerancia_mentor = 7;            % si no mejora en estas rondas, cambia mentor
historial = 10;                   % últimas decisiones que se toman en cuenta para la reputación

% Nota: seguidores + detractores + independientes

%% ---------- Inicialización (no mover) ----------
numero_variables = numel(limites_superiores); 
numero_vecinos = max(2, min(6, floor(numero_variables/3)));   % vecinos por lado en el anillo
rango = limites_superiores - limites_inferiores;
reputacion = 0.5*ones(tamanio_poblacion,historial); % Reputación individual (éxito reciente; EMA de tasa de aciertos)
no_mejora = zeros(tamanio_poblacion,1); % Contadores de no-mejora
vecinos = zeros(tamanio_poblacion,2*numero_vecinos);
convergencia = zeros(maximo_iteraciones,1);
mejorado = false(tamanio_poblacion,1);
ponderacion = (1:historial)/sum(1:historial);
matriz_confianza = ones(size(vecinos)) / (2 * numero_vecinos);
individuos_elegidos = zeros(1,numero_variables);
%% ---------- Población inicial ----------
poblacion = repmat(limites_inferiores, tamanio_poblacion, 1) ...
    + rand(tamanio_poblacion,numero_variables).*repmat(rango, tamanio_poblacion, 1);
costo = zeros(tamanio_poblacion,1);
for i=1:tamanio_poblacion
    costo(i)=funcion_costo(poblacion(i,:)); 
end
[mejor_costo, ibest] = min(costo);
mejor_x = poblacion(ibest,:);

% Memorias sociales
for individuo=1:tamanio_poblacion
    idx = mod((individuo-1)+(-numero_vecinos:numero_vecinos), tamanio_poblacion)+1; 
    idx(idx==individuo)=[];
    vecinos(individuo,:) = idx;
end

%% ---------- Proceso iterativo ----------
for iter = 1:maximo_iteraciones

    % ---- Fase 1: todos proponen una dirección ----
    % Combina: (a) dirección propia aleatoria; (b) sesgo hacia (best_x - X(i,:))

    % ---- Fase 2: cada persona elige a quién seguir (softmax por confianza*reputación) ----
    for individuo=1:tamanio_poblacion
        idx = vecinos(individuo,:); 
        reput = sum(reputacion(idx,:).*repmat(ponderacion,length(idx),1),2)';           % reputación de vecinos
        score = matriz_confianza(individuo,:) .* reput; 
        % incluye la opción “me sigo a mí” como otro "vecino" artificial
        score_self = sum(reputacion(individuo,:) .* ponderacion); 
        idx_aug = [idx(score>score_self), individuo];
        score_aug = [score(score>score_self), score_self];
        % softmax con temperatura
        z = score_aug / temperatura_social;
        % muestreo discreto
        z = z - max(z);
        p = exp(z);
        p = max(p, 1e-12); 
        p = p / sum(p);
        for no_var = 1:numero_variables
            k = find(rand<=cumsum(p),1);
            individuos_elegidos(no_var) = k;
        end
        respuesta_copiada = poblacion(sub2ind(size(poblacion), individuos_elegidos, 1:numero_variables));
    end

    % ---- Fase 3: Stand-up (cada m_standup): líder propone y hay seguidores/retractores ----
    if mod(iter, pasos_reunion) == 0
        [~, ibest] = min(costo); 
        lider = ibest;
        if norm(mejor_x - poblacion(lider,:))>0
            direccion_lider = (mejor_x - poblacion(lider,:)) / norm(mejor_x - poblacion(lider,:));
        else
            direccion_lider = randn(1,numero_variables); % respaldo
            direccion_lider = direccion_lider / norm(direccion_lider);
        end

        % etiquetar seguidores y retractores
        idx_all = randperm(tamanio_poblacion);
        n_seguidores = round(fraccion_seguidores * tamanio_poblacion);
        n_retractores  = round(fraccion_detractores   * tamanio_poblacion);
        lista_seguidores = idx_all(1:n_seguidores);
        lista_retractores  = idx_all(n_seguidores+1:min(n_seguidores+n_retractores, tamanio_poblacion));

        % aplicar ajuste de rumbo
        for seguidor = lista_seguidores
            direccion_elegida(seguidor,:) = (1-sesgo_al_mejor)*direccion_elegida(seguidor,:) + sesgo_al_mejor*direccion_lider;  % se alinean
            direccion_elegida(seguidor,:) = direccion_elegida(seguidor,:) / norm(direccion_elegida(seguidor,:));
        end
        for retractor = lista_retractores
            direccion_elegida(retractor,:) = (1-sesgo_al_mejor)*direccion_elegida(retractor,:) - sesgo_al_mejor*direccion_lider;  % prueban lo contrario
            direccion_elegida(retractor,:) = direccion_elegida(retractor,:) / norm(direccion_elegida(retractor,:));
        end
    end

    % ---- Fase 4: Moverse, evaluar y actualizar confianza/ reputación/ pasos ----
    nueva_poblacion = poblacion; 
    nuevo_costo = costo;

    for individuo=1:tamanio_poblacion
        paso = sigma(individuo,:) .* direccion_elegida(individuo,:);
        candidato = reflect_box(poblacion(individuo,:) + paso, limites_inferiores, limites_superiores);
        costo_candidato = funcion_costo(candidato);

        if costo_candidato < costo(individuo)   % éxito
            nueva_poblacion(individuo,:) = candidato; 
            nuevo_costo(individuo) = costo_candidato; 
            mejorado(individuo) = true;
            sigma(individuo,:) = 1.1 * sigma(individuo,:);
            no_mejora(individuo) = 0;
        else                % fracaso
            sigma(individuo,:) = 0.9 * sigma(individuo,:);
            no_mejora(individuo) = no_mejora(individuo) + 1;
            mejorado(individuo) = false;
        end

        % reputación individual (EMA de éxitos 0/1)
        reputacion(individuo,:) = [reputacion(individuo,2:end),double(mejorado(individuo))];

        % refuerzo de confianza si seguí a un mentor real
        if mentor_elegido(individuo)
            idx = vecinos(individuo,:); 
            idx_mentor = find(idx==mentor(individuo), 1);
            if ~isempty(idx_mentor)
                recompensa = (costo_candidato < costo(individuo)); % recompensa 1 si mejoró
                % introducir un pequeño olvido
                matriz_confianza(individuo,:) = 0.9 * matriz_confianza(individuo,:);
                matriz_confianza(individuo,idx_mentor) = (1-tasa_confianza)*matriz_confianza(individuo,idx_mentor) + tasa_confianza*double(recompensa);
                % renormaliza (mantén pequeña masa a los no-pos)
                matriz_confianza(individuo,:) = matriz_confianza(individuo,:) / sum(matriz_confianza(individuo,:));
            end
        end

        % si llevo mucho sin mejorar cambio de mentor (re-cableo ligero)
        if no_mejora(individuo) >= tolerancia_mentor
            % suelta al mentor actual y añade un atajo aleatorio
            mentor(individuo) = 0;
            no_mejora(individuo) = 0; % reset
        end
    end

    % aplicar cambios
    poblacion = nueva_poblacion; 
    costo = nuevo_costo;

    % actualiza mejor global
    [cmin, imin] = min(costo);
    if cmin < mejor_costo
        mejor_costo = cmin; 
        mejor_x = poblacion(imin,:);
    end

    convergencia(iter) = mejor_costo;
%     plot(X(:,1), X(:,2), 'ok')
%     plot([0],[0], 'xr', 'MarkerSize', 12)
%     hold on
%     xlim([-600,600])
%     ylim([-600,600])
%     pause(0.5)
end
solucion = [mejor_x, mejor_costo];

end 

% ---------- Utilidad: reflexión en caja ----------
function x = reflect_box(x, lo, hi)
below = x < lo; 
if any(below)
    x(below) = 2*lo(below) - x(below);
end
above = x > hi; 
if any(above)
    x(above) = 2*hi(above) - x(above); 
end
x = min(max(x, lo), hi);
end

