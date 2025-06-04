%--------------------------------------------------------------------------
% This script loads per-participant correlation data for “style vs. liking”
% and “category vs. liking” at two different frequencies (20 Hz and 3 Hz).
% It splits “style” correlations into three predefined groups, accumulates
% sums across participants, and computes grand-average timecourses. 
% Finally, it plots the averaged correlations (with statistical markers) 
% for each style group and each category, side by side for 20 Hz and 3 Hz.
%--------------------------------------------------------------------------

clear all;
clc;

%%============================ PARAMETERS ================================

% Number of participants (indexed 18:39 for 20 Hz; 1:39 for 3 Hz)
numparticipants    = 39;
numparticipants_3hz = 39;

% Preallocate accumulators for style (3 groups × 110 timepoints)
sum_style_liking       = zeros(3, 110);
sum_style_liking_3hz   = zeros(3, 110);

% Preallocate accumulators for category (5 categories × 110 timepoints)
sum_category_liking    = zeros(5, 110);
sum_category_liking_3hz = zeros(5, 110);

% Counters for how many files were successfully processed
total_processed       = 0;
total_processed_3hz   = 0;

% Preallocate per-participant storage for “style” group means at 20 Hz
% 22 = maximum possible participants processed (39−18+1 = 22)
style_1    = zeros(22, 110);
style_2    = zeros(22, 110);
style_3    = zeros(22, 110);

% Preallocate per-participant storage for “category” correlations at 20 Hz
cat_1 = [];  % will grow to [N_participants × 110]
cat_2 = [];
cat_3 = [];
cat_4 = [];
cat_5 = [];

% Same preallocations for 3 Hz
style_1_3hz = zeros(34, 110);
style_2_3hz = zeros(34, 110);
style_3_3hz = zeros(34, 110);

cat_1_3hz = [];
cat_2_3hz = [];
cat_3_3hz = [];
cat_4_3hz = [];
cat_5_3hz = [];

% Define which rows of “s” correspond to each style group
group1_indices = [1, 2, 3, 4, 9];   % e.g., “Classical Tendencies” 
group2_indices = [5, 6, 7, 10];     % e.g., “Romantic/Emotive Tendencies”
group3_indices = [8, 11, 12, 13];   % e.g., “Modernist Tendencies”

% Folder where “results_sty_cat_*.mat” files reside
foldername = 'D:\Dropbox\Internship\tasks_may2024\Results\generated_results\sty_cat_result';

%%========================= PROCESS 20 Hz DATA =============================

for i = 18:numparticipants
    % Construct the filename for participant i at 20 Hz
    filename = fullfile(foldername, sprintf('results_sty_cat_20hz_0%d.mat', i));
    
    if exist(filename, 'file') == 2
        data = load(filename);
        % s = 3×110 “style vs. liking” correlations
        % c = 5×110 “category vs. liking” correlations
        corr_style_liking    = data.s;  
        corr_category_liking = data.c;  
        
        % Compute the mean across predefined style groups
        group1 = mean(corr_style_liking(group1_indices, :), 1);  % 1×110
        group2 = mean(corr_style_liking(group2_indices, :), 1);  % 1×110
        group3 = mean(corr_style_liking(group3_indices, :), 1);  % 1×110
        
        % Store each group’s timecourse in a growing array
        p = total_processed + 1;  % current row index
        style_1(p, :) = group1;
        style_2(p, :) = group2;
        style_3(p, :) = group3;
        
        % Append each category’s 1×110 row
        cat_1 = [cat_1; corr_category_liking(1, :)];
        cat_2 = [cat_2; corr_category_liking(2, :)];
        cat_3 = [cat_3; corr_category_liking(3, :)];
        cat_4 = [cat_4; corr_category_liking(4, :)];
        cat_5 = [cat_5; corr_category_liking(5, :)];
        
        % Accumulate sums for grand-average later
        sum_style_liking    = sum_style_liking    + [group1; group2; group3];
        sum_category_liking = sum_category_liking + corr_category_liking;
        
        fprintf('Loaded and processed: %s\n', filename);
        total_processed = total_processed + 1;
    else
        fprintf('File not found: %s\n', filename);
    end
end

%%========================= PROCESS 3 Hz DATA ==============================

for i = 1:numparticipants_3hz
    % Construct the filename for participant i at 3 Hz
    filename = fullfile(foldername, sprintf('results_sty_cat_0%d.mat', i));

    if exist(filename, 'file') == 2
        data = load(filename);
        corr_style_liking    = data.s;  
        corr_category_liking = data.c;  
        
        % Compute mean across the same style groups
        group1 = mean(corr_style_liking(group1_indices, :), 1);
        group2 = mean(corr_style_liking(group2_indices, :), 1);
        group3 = mean(corr_style_liking(group3_indices, :), 1);
        
        % Store into the 3 Hz arrays
        p3 = total_processed_3hz + 1;  % current row index
        style_1_3hz(p3, :) = group1;
        style_2_3hz(p3, :) = group2;
        style_3_3hz(p3, :) = group3;
        
        % Append category rows
        cat_1_3hz = [cat_1_3hz; corr_category_liking(1, :)];
        cat_2_3hz = [cat_2_3hz; corr_category_liking(2, :)];
        cat_3_3hz = [cat_3_3hz; corr_category_liking(3, :)];
        cat_4_3hz = [cat_4_3hz; corr_category_liking(4, :)];
        cat_5_3hz = [cat_5_3hz; corr_category_liking(5, :)];
        
        % Accumulate sums
        sum_style_liking_3hz    = sum_style_liking_3hz    + [group1; group2; group3];
        sum_category_liking_3hz = sum_category_liking_3hz + corr_category_liking;
        
        fprintf('Loaded and processed: %s\n', filename);
        total_processed_3hz = total_processed_3hz + 1;
    else
        fprintf('File not found: %s\n', filename);
    end
end

%%=========================== COMPUTE AVERAGES =============================

% Grand-average for 20 Hz (3 style groups × 110)
average_style    = sum_style_liking    / total_processed;
% Grand-average for 20 Hz (5 categories × 110)
average_cat      = sum_category_liking / total_processed;

% Grand-average for 3 Hz
average_style_3hz = sum_style_liking_3hz    / total_processed_3hz;
average_cat_3hz   = sum_category_liking_3hz / total_processed_3hz;

%%--------------------------- PLOT SETTINGS -------------------------------

% Define style and category labels
style_variables = {'Classical Tendencies', 'Romantic/Emotive Tendencies', 'Modernist Tendencies'};
cat_variables   = {'Scenes', 'Portrait', 'Landscape', 'Still life', 'Toward Abstraction'};

% Generate colormaps
style_colors = lines(numel(style_variables));  % returns num_styles × 3 RGB
cat_colors   = lines(numel(cat_variables));    % returns num_cat × 3 RGB

% Convert to cell arrays of RGB triplets
style_color_cells = mat2cell(style_colors, ones(1, numel(style_variables)), 3);
cat_color_cells   = mat2cell(cat_colors,   ones(1, numel(cat_variables)),   3);

% Create maps from label → RGB
feature_style_colors = containers.Map(style_variables, style_color_cells);
feature_cat_colors   = containers.Map(cat_variables,   cat_color_cells);

% Line widths and styles
lineWidths = [3, 2, 1];  % [highlighted-feature width, other-feature width, marker width]

% Load the universal time vector (1×110) for x-axis
time_data = load('D:\Dropbox\Internship\tasks_may2024\time.mat');
time_vec  = time_data.time;  % assume variable is called “time” inside .mat

%%======================= PLOT: STYLE (Liking) =============================

figure;
t = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% ------------------ Left tile: 20 Hz (50 ms) -------------------
nexttile;
hold on;
legend_handles = [];

% Plot each style’s average timecourse at 20 Hz
for idx = 1:numel(style_variables)
    label    = style_variables{idx};
    rgb_full = feature_style_colors(label);
    
    h = plot(time_vec, average_style(idx, :), ...
        'Color', rgb_full, ...
        'LineStyle', '-', ...
        'LineWidth', lineWidths(2), ...
        'DisplayName', label);
    
    legend_handles = [legend_handles, h];
end

% Shade stimulus-on window (0 to 0.050 s)
p = fill([0, 0.050, 0.050, 0], [-0.007, -0.007, -0.0065, -0.0065], ...
    [0.5, 0.5, 0.5], 'EdgeColor', 'k');
p.LineWidth = 2;

% Dummy patch for “Stimulus on” legend entry
h_stim = patch(NaN, NaN, [0.5, 0.5, 0.5], 'DisplayName', 'Stimulus on', 'LineWidth', 1);
legend_handles = [legend_handles, h_stim];

% Perform one-sample right-tailed t-test at each timepoint per style group
feature_names = {'style_1', 'style_2', 'style_3'};
min_y = min(average_style(:));
max_y = max(average_style(:));
p_adj_matrix = [];

for i = 1:numel(feature_names)
    data_matrix = eval(feature_names{i});  % [N_participants × 110]
    [~, pv]     = ttest(data_matrix, 0, 'tail', 'right');
    [~, ~, ~, adj_pv] = fdr_bh(pv);         % FDR-corrected p-values
    p_adj_matrix = [p_adj_matrix; adj_pv];
    
    % Mark significant timepoints (adj_pv < 0.05 & time > 0)
    sig_points = (adj_pv < 0.05) & (time_vec > 0);
    offset     = (i - 1) * 0.002;  % vertical offset so markers don’t overlap
    yline_val  = max_y + offset + 0.002;
    
    if any(sig_points)
        marker_color = feature_style_colors(style_variables{i});
        plot(time_vec(sig_points), repmat(yline_val, 1, sum(sig_points)), 's', ...
            'MarkerSize', 4, ...
            'MarkerEdgeColor', marker_color, ...
            'MarkerFaceColor', marker_color, ...
            'LineWidth', 1);
    end
end

% Axes decorations
xline(0, '--');
yline(0, '--');
xlim([min(time_vec), max(time_vec)]);
xlabel('Time (s)');
ylabel('Correlation');
title('50 ms Liking for Different Styles (20 Hz)');
legend(legend_handles, 'Location', 'best');
grid off;
hold off;

% ------------------ Right tile: 3 Hz (150 ms) -------------------
nexttile;
hold on;
legend_handles_3hz = [];

for idx = 1:numel(style_variables)
    label    = style_variables{idx};
    rgb_full = feature_style_colors(label);
    
    h = plot(time_vec, average_style_3hz(idx, :), ...
        'Color', rgb_full, ...
        'LineStyle', '-', ...
        'LineWidth', lineWidths(2), ...
        'DisplayName', label);
    
    legend_handles_3hz = [legend_handles_3hz, h];
end

% Shade stimulus-on window (0 to 0.150 s)
p = fill([0, 0.150, 0.150, 0], [-0.007, -0.007, -0.0065, -0.0065], ...
    [0.5, 0.5, 0.5], 'EdgeColor', 'k');
p.LineWidth = 2;

h_stim = patch(NaN, NaN, [0.5, 0.5, 0.5], 'DisplayName', 'Stimulus on', 'LineWidth', 1);
legend_handles_3hz = [legend_handles_3hz, h_stim];

feature_names_3hz = {'style_1_3hz', 'style_2_3hz', 'style_3_3hz'};
min_y3 = min(average_style_3hz(:));
max_y3 = max(average_style_3hz(:));
p_adj_matrix_3hz = [];

for i = 1:numel(feature_names_3hz)
    data_matrix_3hz = eval(feature_names_3hz{i});
    [~, pv3]        = ttest(data_matrix_3hz, 0, 'tail', 'right');
    [~, ~, ~, adj_pv3] = fdr_bh(pv3);
    p_adj_matrix_3hz = [p_adj_matrix_3hz; adj_pv3];
    
    sig_points_3hz = (adj_pv3 < 0.05) & (time_vec > 0);
    offset3        = (i - 1) * 0.002;
    yline_val3     = max_y3 + offset3 + 0.002;
    
    if any(sig_points_3hz)
        marker_color = feature_style_colors(style_variables{i});
        plot(time_vec(sig_points_3hz), repmat(yline_val3, 1, sum(sig_points_3hz)), 's', ...
            'MarkerSize', 4, ...
            'MarkerEdgeColor', marker_color, ...
            'MarkerFaceColor', marker_color, ...
            'LineWidth', 1);
    end
end

xline(0, '--');
yline(0, '--');
xlim([min(time_vec), max(time_vec)]);
xlabel('Time (s)');
ylabel('Correlation');
title('150 ms Liking for Different Styles (3 Hz)');
legend(legend_handles_3hz, 'Location', 'best');
grid off;
hold off;

% Save the style plots
plotpath_style = fullfile(foldername, 'average', 'Style_Liking_correlation.png');
set(gcf, 'Position', [100, 100, 1200, 600]);
exportgraphics(gcf, plotpath_style, 'Resolution', 300);

%%==================== PLOT: CATEGORY (Liking) =============================

figure;
t2 = tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

% ------------------ Left tile: 20 Hz (50 ms) -------------------
nexttile;
hold on;
legend_handles_cat = [];

for idx = 1:numel(cat_variables)
    label    = cat_variables{idx};
    rgb_full = feature_cat_colors(label);
    
    h = plot(time_vec, average_cat(idx, :), ...
        'Color', rgb_full, ...
        'LineStyle', '-', ...
        'LineWidth', lineWidths(2), ...
        'DisplayName', label);
    
    legend_handles_cat = [legend_handles_cat, h];
end

% Shade stimulus-on window (0 to 0.050 s)
p = fill([0, 0.050, 0.050, 0], [-0.007, -0.007, -0.0065, -0.0065], ...
    [0.5, 0.5, 0.5], 'EdgeColor', 'k');
p.LineWidth = 2;

h_stim = patch(NaN, NaN, [0.5, 0.5, 0.5], 'DisplayName', 'Stimulus on', 'LineWidth', 1);
legend_handles_cat = [legend_handles_cat, h_stim];

feature_names_cat = {'cat_1', 'cat_2', 'cat_3', 'cat_4', 'cat_5'};
min_catY = min(average_cat(:));
max_catY = max(average_cat(:));
p_adj_cat = [];

for i = 1:numel(feature_names_cat)
    data_cat = eval(feature_names_cat{i});  % [N_participants × 110]
    [~, pv_cat] = ttest(data_cat, 0, 'tail', 'right');
    [~, ~, ~, adj_pv_cat] = fdr_bh(pv_cat);
    p_adj_cat = [p_adj_cat; adj_pv_cat];
    
    sig_points_cat = (adj_pv_cat < 0.05) & (time_vec > 0);
    offset_cat     = (i - 1) * 0.002;
    yline_cat      = max_catY + offset_cat + 0.002;
    
    if any(sig_points_cat)
        marker_color = feature_cat_colors(cat_variables{i});
        plot(time_vec(sig_points_cat), repmat(yline_cat, 1, sum(sig_points_cat)), 's', ...
            'MarkerSize', 4, ...
            'MarkerEdgeColor', marker_color, ...
            'MarkerFaceColor', marker_color, ...
            'LineWidth', 1);
    end
end

xline(0, '--');
yline(0, '--');
xlim([min(time_vec), max(time_vec)]);
xlabel('Time (s)');
ylabel('Correlation');
title('50 ms Liking for Different Categories (20 Hz)');
legend(legend_handles_cat, 'Location', 'best');
grid off;
hold off;

% ------------------ Right tile: 3 Hz (150 ms) -------------------
nexttile;
hold on;
legend_handles_cat_3hz = [];

for idx = 1:numel(cat_variables)
    label    = cat_variables{idx};
    rgb_full = feature_cat_colors(label);
    
    h = plot(time_vec, average_cat_3hz(idx, :), ...
        'Color', rgb_full, ...
        'LineStyle', '-', ...
        'LineWidth', lineWidths(2), ...
        'DisplayName', label);
    
    legend_handles_cat_3hz = [legend_handles_cat_3hz, h];
end

% Shade stimulus-on window (0 to 0.150 s)
p = fill([0, 0.150, 0.150, 0], [-0.007, -0.007, -0.0065, -0.0065], ...
    [0.5, 0.5, 0.5], 'EdgeColor', 'k');
p.LineWidth = 2;

h_stim = patch(NaN, NaN, [0.5, 0.5, 0.5], 'DisplayName', 'Stimulus on', 'LineWidth', 1);
legend_handles_cat_3hz = [legend_handles_cat_3hz, h_stim];

feature_names_cat_3hz = {'cat_1_3hz', 'cat_2_3hz', 'cat_3_3hz', 'cat_4_3hz', 'cat_5_3hz'};
min_catY3 = min(average_cat_3hz(:));
max_catY3 = max(average_cat_3hz(:));
p_adj_cat_3hz = [];

for i = 1:numel(feature_names_cat_3hz)
    data_cat_3hz = eval(feature_names_cat_3hz{i});
    [~, pv_cat_3hz] = ttest(data_cat_3hz, 0, 'tail', 'right');
    [~, ~, ~, adj_pv_cat_3hz] = fdr_bh(pv_cat_3hz);
    p_adj_cat_3hz = [p_adj_cat_3hz; adj_pv_cat_3hz];
    
    sig_points_cat_3hz = (adj_pv_cat_3hz < 0.05) & (time_vec > 0);
    offset_cat_3hz     = (i - 1) * 0.002;
    yline_cat_3hz      = max_catY3 + offset_cat_3hz + 0.002;
    
    if any(sig_points_cat_3hz)
        marker_color = feature_cat_colors(cat_variables{i});
        plot(time_vec(sig_points_cat_3hz), repmat(yline_cat_3hz, 1, sum(sig_points_cat_3hz)), 's', ...
            'MarkerSize', 4, ...
            'MarkerEdgeColor', marker_color, ...
            'MarkerFaceColor', marker_color, ...
            'LineWidth', 1);
    end
end

xline(0, '--');
yline(0, '--');
xlim([min(time_vec), max(time_vec)]);
xlabel('Time (s)');
ylabel('Correlation');
title('150 ms Liking for Different Categories (3 Hz)');
legend(legend_handles_cat_3hz, 'Location', 'best');
grid off;
hold off;

% Save the category plots
plotpath_cat = fullfile(foldername, 'average', 'Category_Liking_correlation.png');
set(gcf, 'Position', [100, 100, 1200, 600]);
exportgraphics(gcf, plotpath_cat, 'Resolution', 300);