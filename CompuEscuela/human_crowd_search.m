function [solucion, convergencia] = human_crowd_search(funcion_costo, limites_inferiores, limites_superiores, tamanio_poblacion, maximo_iteraciones)
% BUSQUEDA_GRUPAL - Human Crowd Search (HCS)
% B�squeda inspirada en grupos humanos:
% - Red social ligera (anillo + atajos)
% - Consejos entre pares (direcciones compartidas)
% - Confianza/reputaci�n que se actualiza por refuerzo
% - Reuniones peri�dicas ("stand-ups") con l�der y disenso controlado
% - Satisficing: cambio de mentor cuando no hay progreso
% - Adaptaci�n de paso por individuo y reflexi�n en l�mites
%
% Entradas: como tu versi�n original.
% Salidas:
%   solucion      : [1 x (d+1)] = [mejor_x, mejor_costo]
%   convergencia  : [maximo_iteraciones x 1] costo m�nimo global por iter

%% ---------- Par�metros de control ----------
fraccion_seguidores = 0.5;           % fracci�n que sigue al l�der en stand-up
fraccion_detractores   = 0.2;           % fracci�n que prueba rumbo opuesto
pasos_reunion = 13;                    % cada cu�ntas iteraciones hay "reuni�n"
temperatura_social = 0.7;                % temperatura para softmax de consejos
tasa_confianza = 0.15;                 % tasa de aprendizaje de confianza
tolerancia_mentor = 7;            % si no mejora en estas rondas, cambia mentor
sesgo_al_mejor = 0.4;             % 40% sesgo al mejor global, 60% exploraci�n propia
historial = 10;                   % �ltimas decisiones que se toman en cuenta para la reputaci�n
sesgo_mentor = 0.3;               % es la ventaja que se le da al mentor de ser elegido

% Nota: seguidores + detractores + independientes

%% ---------- Inicializaci�n (no mover) ----------
numero_variables = numel(limites_superiores); 
numero_vecinos = max(2, min(6, floor(numero_variables/3)));   % vecinos por lado en el anillo
rango = limites_superiores - limites_inferiores;
reputacion = 0.5*ones(tamanio_poblacion,historial); % Reputaci�n individual (�xito reciente; EMA de tasa de aciertos)
no_mejora = zeros(tamanio_poblacion,1); % Contadores de no-mejora
vecinos = zeros(tamanio_poblacion,2*numero_vecinos);
mentor = zeros(tamanio_poblacion,1); % Mentor preferido por cada individuo (arranca sin mentor: 0 => s�guete a ti)
sigma0 = rango / 20;
convergencia = zeros(maximo_iteraciones,1);
direccion_elegida = zeros(tamanio_poblacion,numero_variables);
mentor_elegido = false(tamanio_poblacion,1);
mejorado = false(tamanio_poblacion,1);
ponderacion = (1:historial)/sum(1:historial);
matriz_confianza = ones(size(vecinos)) / (2 * numero_vecinos);

%% ---------- Poblaci�n inicial ----------
poblacion = repmat(limites_inferiores, tamanio_poblacion, 1) ...
    + rand(tamanio_poblacion,numero_variables).*repmat(rango, tamanio_poblacion, 1);
sigma = repmat(sigma0, tamanio_poblacion, 1);                  % paso por componente
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

    % ---- Fase 1: todos proponen una direcci�n ----
    % Combina: (a) direcci�n propia aleatoria; (b) sesgo hacia (best_x - X(i,:))
    direcciones = zeros(tamanio_poblacion,numero_variables);
    for individuo=1:tamanio_poblacion
        vector_unitario = randn(1,numero_variables); 
        vector_unitario = vector_unitario / norm(vector_unitario);
        sesgo = mejor_x - poblacion(individuo,:);
        tamanio_sesgo = norm(sesgo);
        if tamanio_sesgo > 0
            sesgo = sesgo / tamanio_sesgo;
            direccion_propuesta = sesgo_al_mejor * sesgo + (1-sesgo_al_mejor) * vector_unitario;
        else
            direccion_propuesta = vector_unitario;
        end
        direcciones(individuo,:) = direccion_propuesta / norm(direccion_propuesta);
    end

    % ---- Fase 2: cada persona elige a qui�n seguir (softmax por confianza*reputaci�n) ----
    for individuo=1:tamanio_poblacion
        idx = vecinos(individuo,:); 
        reput = sum(reputacion(idx,:).*repmat(ponderacion,length(idx),1),2)';           % reputaci�n de vecinos
        score = matriz_confianza(individuo,:) .* reput; 
        % incluye la opci�n �me sigo a m� como otro "vecino" artificial
        score_self = sum(reputacion(individuo,:) .* ponderacion); 
        idx_aug = [idx, tamanio_poblacion+individuo];   % n+i representar� "self"
        score_aug = [score, score_self];
        % softmax con temperatura
        z = score_aug / temperatura_social;
        % si tengo mentor fijado y conf�o en �l, aumento su probabilidad
        if mentor(individuo)>0
            k = find(idx_aug==mentor(individuo), 1);
            % Reajustar esta parte, est� horrible
            if ~isempty(k)
                z(k) = z(k) + sesgo_mentor;
            end
        end
        % muestreo discreto
        z = z - max(z);
        p = exp(z);
        p = max(p, 1e-12); 
        p = p / sum(p);
        k = find(rand<=cumsum(p),1);
        individuo_elegido = idx_aug(k);
        if individuo_elegido <= tamanio_poblacion
            direccion_elegida(individuo,:) = direcciones(individuo_elegido,:); % sigo consejo del vecino
            mentor_elegido(individuo) = true;
            mentor(individuo) = individuo_elegido;               % actualizo mentor preferido
        else
            direccion_elegida(individuo,:) = direcciones(individuo,:);      % sigo mi propio consejo
            mentor_elegido(individuo) = false;
            mentor(individuo) = 0;
        end
    end

    % ---- Fase 3: Stand-up (cada m_standup): l�der propone y hay seguidores/retractores ----
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

    % ---- Fase 4: Moverse, evaluar y actualizar confianza/ reputaci�n/ pasos ----
    nueva_poblacion = poblacion; 
    nuevo_costo = costo;

    for individuo=1:tamanio_poblacion
        paso = sigma(individuo,:) .* direccion_elegida(individuo,:);
        candidato = reflect_box(poblacion(individuo,:) + paso, limites_inferiores, limites_superiores);
        costo_candidato = funcion_costo(candidato);

        if costo_candidato < costo(individuo)   % �xito
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

        % reputaci�n individual (EMA de �xitos 0/1)
        reputacion(individuo,:) = [reputacion(individuo,2:end),double(mejorado(individuo))];

        % refuerzo de confianza si segu� a un mentor real
        if mentor_elegido(individuo)
            idx = vecinos(individuo,:); 
            idx_mentor = find(idx==mentor(individuo), 1);
            if ~isempty(idx_mentor)
                recompensa = (costo_candidato < costo(individuo)); % recompensa 1 si mejor�
                % introducir un peque�o olvido
                matriz_confianza(individuo,:) = 0.9 * matriz_confianza(individuo,:);
                matriz_confianza(individuo,idx_mentor) = (1-tasa_confianza)*matriz_confianza(individuo,idx_mentor) + tasa_confianza*double(recompensa);
                % renormaliza (mant�n peque�a masa a los no-pos)
                matriz_confianza(individuo,:) = matriz_confianza(individuo,:) / sum(matriz_confianza(individuo,:));
            end
        end

        % si llevo mucho sin mejorar cambio de mentor (re-cableo ligero)
        if no_mejora(individuo) >= tolerancia_mentor
            % suelta al mentor actual y a�ade un atajo aleatorio
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

% ---------- Utilidad: reflexi�n en caja ----------
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
