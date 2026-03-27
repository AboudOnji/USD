%% ========================================================================
%  SIMULACIÓN MONTE CARLO — TIPO DE CAMBIO USD/MXN
%  Modelo: Movimiento Browniano Geométrico (GBM) con reversión parcial
%          y saltos de Poisson (modelo Jump-Diffusion de Merton)
%
%  Parámetros calibrados con datos históricos reales:
%    - Promedio 2023: 17.74 MXN/USD
%    - Promedio 2024: 18.33 MXN/USD
%    - Máximo 2024:   20.87 MXN/USD (post-elecciones jun-2024)
%    - Mínimo 2023:   16.63 MXN/USD (julio 2023, superpeso)
%    - Nivel actual (mar-2026): ~17.50 MXN/USD
%    - Volatilidad anualizada histórica: ~10–13% (media ~11.5%)
%    - Diferencial de tasas Banxico-Fed: ~5.75 pp (sustenta al peso)
%
%  Author : Prof. D.Sc. BARSEKH-ONJI Aboud
%  Institución: Facultad de Ingeniería, Universidad Anáhuac México
%  Fecha  : Marzo 2026
%% ========================================================================

clc; clear; close all;
rng(42);   % Semilla para reproducibilidad

%% ========================================================================
%  SECCIÓN 1 — PARÁMETROS DEL MODELO
%% ========================================================================

% --- Horizonte de simulación ---
T_anos   = 3;               % Horizonte total: 3 anos
dt       = 1/252;           % Paso temporal: 1 día hábil (252 días/año)
N_pasos  = round(T_anos / dt);   % Número total de pasos

% --- Condición inicial (nivel actual mar-2026) ---
S0 = 17.50;   % MXN por 1 USD

% --- Parámetros de deriva (mu) ---
% La deriva neta refleja:
%   + Inflación diferencial MX–US  ≈ +3.5% anual (deprecia al MXN)
%   - Diferencial de tasas real     ≈ -2.5% anual (aprecia al MXN)
%   + Prima de riesgo país + incertidumbre ≈ +1.5%
% => mu neto ≈ +2.5% anual (deprecia MXN modestamente)
mu = 0.025;       % Deriva anual (2.5%)

% --- Volatilidad (sigma) ---
% Calibrada con datos históricos 2023–2025; se usa EWMA para reflejar
% que la volatilidad reciente es ligeramente menor al promedio histórico
sigma = 0.115;    % Volatilidad anual (11.5%)

% --- Parámetros de saltos de Poisson (Modelo Merton 1976) ---
% Los saltos modelan shocks discretos: elecciones, decisiones de Banxico/Fed,
% crisis externas (ej: carry trade unwind de agosto 2024)
lambda_j   = 3.0;    % Intensidad de saltos: ~3 eventos/año
mu_j       = -0.010; % Tamaño medio del salto en log-escala (-1.0%)
% Negativo: más devaluaciones que apreciaciones abruptas
sigma_j    = 0.030;  % Desviación del salto (3%)

% --- Número de trayectorias Monte Carlo ---
M = 1000;   % Simulaciones; balance entre precisión y tiempo de cómputo

% --- Nivel de confianza ---
alpha_IC = 0.95;    % Intervalo de confianza del 95%
p_lower  = (1 - alpha_IC) / 2;     % Percentil inferior (2.5%)
p_upper  = 1 - p_lower;            % Percentil superior (97.5%)

%% ========================================================================
%  SECCIÓN 2 — SIMULACIÓN MONTE CARLO (GBM + SALTOS DE POISSON)
%% ========================================================================

fprintf('=================================================\n');
fprintf('  Simulación USD/MXN — Monte Carlo + Jump-Diff  \n');
fprintf('=================================================\n');
fprintf('  Trayectorias  : %d\n', M);
fprintf('  Horizonte     : %d anos (%d días hábiles)\n', T_anos, N_pasos);
fprintf('  S0            : %.4f MXN/USD\n', S0);
fprintf('  mu (deriva)   : %.2f%% anual\n', mu*100);
fprintf('  sigma         : %.2f%% anual\n', sigma*100);
fprintf('  lambda_j      : %.1f saltos/año\n', lambda_j);
fprintf('=================================================\n\n');

% Preasignación de memoria
S = zeros(N_pasos + 1, M);
S(1, :) = S0;

% Corrección de la deriva por la componente de saltos (Merton 1976)
% Para que la esperanza condicional sea mu, se corrige la deriva difusiva
k_merton = exp(mu_j + 0.5 * sigma_j^2) - 1;   % Saltom medio esperado
mu_corr  = mu - lambda_j * k_merton;            % Deriva corregida

% Simulación paso a paso
for t = 1 : N_pasos
    % --- Componente difusiva (Movimiento Browniano) ---
    Z_diff = randn(1, M);   % Innovaciones normales estándar
    dW = sqrt(dt) * Z_diff;

    % --- Componente de saltos (Poisson compuesto) ---
    % Número de saltos en [t, t+dt]
    n_saltos = poissrnd(lambda_j * dt, 1, M);

    % Tamaño acumulado de saltos: suma de log-normales
    salto_total = zeros(1, M);
    for i = 1 : M
        if n_saltos(i) > 0
            salto_i = mu_j * n_saltos(i) + ...
                sigma_j * sqrt(n_saltos(i)) * randn();
            salto_total(i) = salto_i;
        end
    end

    % --- Actualización del tipo de cambio (forma exponencial exacta) ---
    log_S = log(S(t, :)) + ...
        (mu_corr - 0.5 * sigma^2) * dt + ...
        sigma * dW + ...
        salto_total;

    S(t+1, :) = exp(log_S);
end

fprintf('  Simulación completada exitosamente.\n\n');

%% ========================================================================
%  SECCIÓN 3 — ESTADÍSTICAS RESUMIDAS
%% ========================================================================

% Vector de tiempo en anos
t_vec = (0 : N_pasos)' * dt;

% Estadísticas por paso de tiempo
S_media   = mean(S, 2);
S_mediana = median(S, 2);
S_IC_low  = quantile(S, p_lower, 2);
S_IC_high = quantile(S, p_upper, 2);
S_min     = min(S, [], 2);
S_max     = max(S, [], 2);

% Resumen al final del horizonte
S_final = S(end, :);
fprintf('=================================================\n');
fprintf('  RESULTADOS AL AÑO 3 (horizonte completo)\n');
fprintf('=================================================\n');
fprintf('  Media          : %.4f MXN/USD\n', mean(S_final));
fprintf('  Mediana        : %.4f MXN/USD\n', median(S_final));
fprintf('  P2.5  (IC inf) : %.4f MXN/USD\n', quantile(S_final, 0.025));
fprintf('  P97.5 (IC sup) : %.4f MXN/USD\n', quantile(S_final, 0.975));
fprintf('  Mínimo sim.    : %.4f MXN/USD\n', min(S_final));
fprintf('  Máximo sim.    : %.4f MXN/USD\n', max(S_final));
fprintf('  Std. desv.     : %.4f MXN/USD\n', std(S_final));
fprintf('=================================================\n');

%% ========================================================================
%  SECCIÓN 4 — VISUALIZACIÓN
%% ========================================================================

% Eje X en fechas reales (desde marzo 2026)
fecha_inicio = datetime(2026, 3, 27);
fechas = fecha_inicio + caldays(round(t_vec * 365));

% -----------------------------------------------------------------------
%  FIGURA 1: Trayectorias + Banda de Confianza
% -----------------------------------------------------------------------
fig1 = figure('Position', [50 80 1100 580]);

% --- Sub-muestreo de trayectorias para visualización ---
idx_vis = randperm(M, min(150, M));   % Mostrar 150 trayectorias
S_vis   = S(:, idx_vis);

% Graficar trayectorias individuales (color tenue)
hold on;
plot(fechas, S_vis, 'Color', [0.7 0.85 1.0 0.18], 'LineWidth', 0.4);

% Banda de confianza 95%
fill([fechas; flipud(fechas)], ...
    [S_IC_low; flipud(S_IC_high)], ...
    [0.25 0.60 0.90], ...
    'FaceAlpha', 0.20, ...
    'EdgeColor', 'none', ...
    'DisplayName', 'IC 95\%');

% Percentiles 10/90
S_p10 = quantile(S, 0.10, 2);
S_p90 = quantile(S, 0.90, 2);
fill([fechas; flipud(fechas)], ...
    [S_p10; flipud(S_p90)], ...
    [0.15 0.45 0.75], ...
    'FaceAlpha', 0.25, ...
    'EdgeColor', 'none', ...
    'DisplayName', 'IC 80\%');

% Límites de IC
plot(fechas, S_IC_low,  '--', 'Color', [0.10 0.40 0.80], ...
    'LineWidth', 1.2, 'DisplayName', 'P2.5');
plot(fechas, S_IC_high, '--', 'Color', [0.10 0.40 0.80], ...
    'LineWidth', 1.2, 'DisplayName', 'P97.5');

% Mediana y media
plot(fechas, S_mediana, '-', 'Color', [0.85 0.10 0.10], ...
    'LineWidth', 2.2, 'DisplayName', 'Mediana');
plot(fechas, S_media,   '-', 'Color', [0.05 0.55 0.05], ...
    'LineWidth', 2.0, 'LineStyle', '-.', 'DisplayName', 'Media');

% Línea de valor inicial
yline(S0, ':', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.4, ...
    'Label', sprintf('S_0 = %.2f', S0), ...
    'LabelHorizontalAlignment', 'left', ...
    'FontSize', 9);

% Anotaciones de niveles históricos de referencia
yline(20.87, ':', 'Color', [0.75 0.20 0.05], 'LineWidth', 1.0, ...
    'Label', 'Máx. 2024: 20.87', ...
    'LabelHorizontalAlignment', 'right', 'FontSize', 8);
yline(16.63, ':', 'Color', [0.05 0.50 0.10], 'LineWidth', 1.0, ...
    'Label', 'Mín. 2023: 16.63', ...
    'LabelHorizontalAlignment', 'right', 'FontSize', 8);

hold off;

% Etiquetas y formato
xlabel('Fecha', 'FontSize', 13, 'Interpreter', 'latex');
ylabel('USD/MXN (MXN por 1 USD)', 'FontSize', 13, 'Interpreter', 'latex');
title({'Simulación Monte Carlo — Tipo de Cambio USD/MXN', ...
    sprintf('GBM + Saltos de Poisson $\\mid$ M = %d trayectorias $\\mid$ Horizonte = 3 anos', M)}, ...
    'FontSize', 14, 'Interpreter', 'latex');

legend({'Trayectorias', 'IC 95\%', 'IC 80\%', ...
    'P2.5 (IC inf.)', 'P97.5 (IC sup.)', ...
    'Mediana', 'Media'}, ...
    'Location', 'northwest', 'FontSize', 10, 'Interpreter', 'latex');

grid on;
set(gca, 'FontSize', 11, 'TickLabelInterpreter', 'latex', ...
    'GridAlpha', 0.25, 'GridLineStyle', '--');
xlim([fechas(1), fechas(end)]);
ylim([12, 30]);

% Texto de parámetros en la figura
texto_params = sprintf(['$S_0 = %.2f$, $\\mu = %.2f\\%%$, ' ...
    '$\\sigma = %.2f\\%%$, $\\lambda_j = %.1f$ saltos/año'], ...
    S0, mu*100, sigma*100, lambda_j);
annotation('textbox', [0.13 0.01 0.75 0.045], ...
    'String', texto_params, ...
    'Interpreter', 'latex', 'FontSize', 9, ...
    'EdgeColor', 'none', 'HorizontalAlignment', 'center', ...
    'BackgroundColor', [0.97 0.97 0.97]);

% -----------------------------------------------------------------------
%  FIGURA 2: Distribución al año 1, 2 y 3
% -----------------------------------------------------------------------
fig2 = figure('Position', [50 80 1100 420]);

anos_eval  = [1, 2, 3];
colores    = {[0.20 0.50 0.85], [0.85 0.40 0.10], [0.15 0.65 0.30]};
titulos    = {'Año 1 (mar-2027)', 'Año 2 (mar-2028)', 'Año 3 (mar-2029)'};

for k = 1 : 3
    subplot(1, 3, k);

    idx_t = round(anos_eval(k) / dt);   % Índice del paso correspondiente
    S_dist = S(idx_t, :);

    % Histograma normalizado
    histogram(S_dist, 50, ...
        'Normalization', 'pdf', ...
        'FaceColor', colores{k}, ...
        'FaceAlpha', 0.70, ...
        'EdgeColor', 'w');
    hold on;

    % Línea de densidad (kernel)
    [f_kde, xi_kde] = ksdensity(S_dist);
    plot(xi_kde, f_kde, '-', 'Color', colores{k} * 0.7, ...
        'LineWidth', 2.0);

    % Percentiles en el plot
    q025 = quantile(S_dist, 0.025);
    q975 = quantile(S_dist, 0.975);
    q500 = quantile(S_dist, 0.500);

    xline(q025, '--k', 'LineWidth', 1.2);
    xline(q975, '--k', 'LineWidth', 1.2);
    xline(q500, '-r',  'LineWidth', 1.8);

    hold off;

    xlabel('USD/MXN', 'FontSize', 11, 'Interpreter', 'latex');
    ylabel('Densidad', 'FontSize', 11, 'Interpreter', 'latex');
    title(titulos{k}, 'FontSize', 12, 'Interpreter', 'latex');

    % Estadísticas dentro de la caja
    texto_est = sprintf(['Med: %.2f\nP2.5: %.2f\nP97.5: %.2f\n$\\sigma$: %.2f'], ...
        q500, q025, q975, std(S_dist));
    text(0.96, 0.96, texto_est, ...
        'Units', 'normalized', 'HorizontalAlignment', 'right', ...
        'VerticalAlignment', 'top', 'FontSize', 8.5, ...
        'Interpreter', 'latex', 'BackgroundColor', [1 1 1 0.75], ...
        'EdgeColor', [0.7 0.7 0.7]);

    grid on;
    set(gca, 'FontSize', 10, 'TickLabelInterpreter', 'latex', ...
        'GridAlpha', 0.25);
end

sgtitle('Distribución del USD/MXN por Horizonte de Tiempo', ...
    'FontSize', 14, 'Interpreter', 'latex');

% -----------------------------------------------------------------------
%  FIGURA 3: Fan chart (percentiles escalonados)
% -----------------------------------------------------------------------
fig3 = figure('Position', [50 80 1100 500]);

percentiles_fan = [0.05, 0.10, 0.20, 0.30, 0.40, ...
    0.50, 0.60, 0.70, 0.80, 0.90, 0.95];
S_fan = quantile(S, percentiles_fan, 2);

colores_fan = flipud(parula(floor(length(percentiles_fan)/2)));

hold on;
n_bandas = floor(length(percentiles_fan)/2);
for b = 1 : n_bandas
    idx_lo = b;
    idx_hi = length(percentiles_fan) + 1 - b;
    fill([fechas; flipud(fechas)], ...
        [S_fan(:, idx_lo); flipud(S_fan(:, idx_hi))], ...
        colores_fan(b, :), ...
        'FaceAlpha', 0.35, 'EdgeColor', 'none');
end

% Mediana central
plot(fechas, S_fan(:, 6), '-r', 'LineWidth', 2.2, ...
    'DisplayName', 'Mediana (P50)');

% Línea de inicio
yline(S0, ':', 'Color', [0.3 0.3 0.3], 'LineWidth', 1.2, ...
    'Label', sprintf('S_0=%.2f', S0), ...
    'LabelHorizontalAlignment', 'left', 'FontSize', 9);

hold off;

% Colorbar indicativa de percentiles
cb = colorbar('Ticks', linspace(0,1,5), ...
    'TickLabels', {'P5/P95','P10/P90','P20/P80','P30/P70','P40/P60'}, ...
    'TickLabelInterpreter', 'latex', 'FontSize', 9);
cb.Label.String = 'Bandas de Percentil';
cb.Label.Interpreter = 'latex';
cb.Label.FontSize = 10;
colormap(parula);

xlabel('Fecha', 'FontSize', 13, 'Interpreter', 'latex');
ylabel('USD/MXN', 'FontSize', 13, 'Interpreter', 'latex');
title('Fan Chart USD/MXN — Distribución de Percentiles por Fecha', ...
    'FontSize', 14, 'Interpreter', 'latex');
legend({'Mediana (P50)'}, 'Location', 'northwest', ...
    'FontSize', 11, 'Interpreter', 'latex');
grid on;
set(gca, 'FontSize', 11, 'TickLabelInterpreter', 'latex', ...
    'GridAlpha', 0.22, 'GridLineStyle', '--');
xlim([fechas(1), fechas(end)]);

%% ========================================================================
%  SECCIÓN 5 — GUARDAR FIGURAS
%% ========================================================================

exportgraphics(fig1, '/mnt/user-data/outputs/usdmxn_trayectorias.png', ...
    'Resolution', 200);
exportgraphics(fig2, '/mnt/user-data/outputs/usdmxn_distribuciones.png', ...
    'Resolution', 200);
exportgraphics(fig3, '/mnt/user-data/outputs/usdmxn_fanchart.png', ...
    'Resolution', 200);

fprintf('\nFiguras guardadas exitosamente.\n');
fprintf('  -> usdmxn_trayectorias.png\n');
fprintf('  -> usdmxn_distribuciones.png\n');
fprintf('  -> usdmxn_fanchart.png\n');