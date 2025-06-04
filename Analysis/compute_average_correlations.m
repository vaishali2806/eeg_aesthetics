%==========================================================================
% Script to compute average simple and partial correlations across participants
% for both 20 Hz and 3 Hz experiments, segregate by before/after conditions 
% for 20 Hz, and perform statistical tests (t-tests with FDR correction) on 
% these averages.Generates plots of partial and simple correlations over time 
% and saves figures.
%==========================================================================

clear all;
clc;

%%---------------------------- Initialization -----------------------------
% Number of participants for 20 Hz and 3 Hz experiments
numparticipants    = 39;   % for 20 Hz files indexed 18:39
numparticipants_3hz = 39;  % for 3 Hz files indexed 1:39

% Preallocate sums for simple and partial correlations (7 features × 110 timepoints)
sum_simple_corr         = zeros(7, 110);
sum_partial_corr        = zeros(7, 110);
sum_simple_corr_before  = zeros(7, 110);
sum_partial_corr_before = zeros(7, 110);
sum_simple_corr_after   = zeros(7, 110);
sum_partial_corr_after  = zeros(7, 110);

sum_simple_corr_3hz  = zeros(7, 110);
sum_partial_corr_3hz = zeros(7, 110);

% Counters for how many files were successfully loaded/processed
total_processed       = 0;
before_processed      = 0;
after_processed       = 0;
total_processed_3hz   = 0;

% Participant indices for before/after conditions in the 20 Hz experiment
before_20hz = [18, 19, 20, 21, 26, 29, 30, 32, 35, 37, 39];
after_20hz  = [22, 23, 23, 25, 27, 28, 31, 33, 34, 36, 38];

foldername = "D:\Dropbox\Internship\tasks_may2024\Results\generated_results";

%%----------------- Load and Sum Correlations for 20 Hz ------------------
for i = 18:numparticipants
    filename = fullfile(foldername, sprintf('results_20hz_0%d.mat', i));
    
    if exist(filename, "file") == 2
        % Load the .mat file which contains:
        %   s = simple correlation matrix (7 × 110)
        %   p = partial correlation matrix (7 × 110)
        data = load(filename);
        corr_behav_eeg_final = data.s;
        partial_corr_behav_eeg = data.p;

        % Accumulate sums across all participants
        sum_simple_corr  = sum_simple_corr  + corr_behav_eeg_final;
        sum_partial_corr = sum_partial_corr + partial_corr_behav_eeg;
        
        % Accumulate sums separately for before and after conditions
        if ismember(i, before_20hz)
            sum_simple_corr_before  = sum_simple_corr_before  + corr_behav_eeg_final;
            sum_partial_corr_before = sum_partial_corr_before + partial_corr_behav_eeg;
            before_processed = before_processed + 1;
            
        elseif ismember(i, after_20hz)
            sum_simple_corr_after  = sum_simple_corr_after  + corr_behav_eeg_final;
            sum_partial_corr_after = sum_partial_corr_after + partial_corr_behav_eeg;
            after_processed = after_processed + 1;
        end
        
        fprintf('Loaded and processed: %s\n', filename);
        total_processed = total_processed + 1;
    else
        fprintf('File not found: %s\n', filename);
    end
end

%%----------------- Load and Sum Correlations for 3 Hz -------------------
for i = 1:numparticipants_3hz
    filename_3hz = fullfile(foldername, sprintf('results_0%d.mat', i));
    
    if exist(filename_3hz, "file") == 2
        data_3hz = load(filename_3hz);
        corr_behav_eeg_final_3hz = data_3hz.s;
        partial_corr_behav_eeg_3hz = data_3hz.p;
        
        sum_simple_corr_3hz  = sum_simple_corr_3hz  + corr_behav_eeg_final_3hz;
        sum_partial_corr_3hz = sum_partial_corr_3hz + partial_corr_behav_eeg_3hz;
        
        fprintf('Loaded and processed: %s\n', filename_3hz);
        total_processed_3hz = total_processed_3hz + 1;
    else
        fprintf('File not found: %s\n', filename_3hz);
    end
end

%%-------------------------- Compute Averages -----------------------------
% Overall averages for 20 Hz
average_simple            = sum_simple_corr  / total_processed;
average_partial           = sum_partial_corr / total_processed;
average_simple_before     = sum_simple_corr_before  / before_processed;
average_partial_before    = sum_partial_corr_before / before_processed;
average_simple_after      = sum_simple_corr_after  / after_processed;
average_partial_after     = sum_partial_corr_after / after_processed;

% Averages for 3 Hz
average_simple_3hz  = sum_simple_corr_3hz  / total_processed_3hz;
average_partial_3hz = sum_partial_corr_3hz / total_processed_3hz;

%%--------------- Collect Individual Participant Data --------------------
% Initialize matrices to store per-participant correlation timecourses
% Each row corresponds to a participant; each column to a timepoint (110)
s_liking    = [];
s_valence   = [];
s_arousal   = [];
s_complexity = [];
s_familiarity = [];
s_style     = [];
s_category  = [];

p_liking    = [];
p_valence   = [];
p_arousal   = [];
p_complexity = [];
p_familiarity = [];
p_style     = [];
p_category  = [];

s_liking_3hz    = [];
s_valence_3hz   = [];
s_arousal_3hz   = [];
s_complexity_3hz = [];
s_familiarity_3hz = [];
s_style_3hz     = [];
s_category_3hz  = [];

p_liking_3hz    = [];
p_valence_3hz   = [];
p_arousal_3hz   = [];
p_complexity_3hz = [];
p_familiarity_3hz = [];
p_style_3hz     = [];
p_category_3hz  = [];

% Loop again to build per-participant feature matrices for 20 Hz
for participant = 18:numparticipants
    filename = fullfile(foldername, sprintf('results_20hz_0%d.mat', participant));
    
    if exist(filename, "file") == 2
        data = load(filename);
        simple_corr  = data.s;  % 7 × 110 matrix
        partial_corr = data.p;  % 7 × 110 matrix
        
        % Append each feature row to the corresponding matrix
        s_liking      = [s_liking;      simple_corr(1, :)];
        s_valence     = [s_valence;     simple_corr(2, :)];
        s_arousal     = [s_arousal;     simple_corr(3, :)];
        s_complexity  = [s_complexity;  simple_corr(4, :)];
        s_familiarity = [s_familiarity; simple_corr(5, :)];
        s_style       = [s_style;       simple_corr(6, :)];
        s_category    = [s_category;    simple_corr(7, :)];
        
        p_liking      = [p_liking;      partial_corr(1, :)];
        p_valence     = [p_valence;     partial_corr(2, :)];
        p_arousal     = [p_arousal;     partial_corr(3, :)];
        p_complexity  = [p_complexity;  partial_corr(4, :)];
        p_familiarity = [p_familiarity; partial_corr(5, :)];
        p_style       = [p_style;       partial_corr(6, :)];
        p_category    = [p_category;    partial_corr(7, :)];
    else
        fprintf('File not found: %s\n', filename);
    end
end

% Loop again to build per-participant feature matrices for 3 Hz
for participant_3hz = 1:numparticipants_3hz
    filename_3hz = fullfile(foldername, sprintf('results_0%d.mat', participant_3hz));
    
    if exist(filename_3hz, "file") == 2
        data_3hz = load(filename_3hz);
        simple_corr_3hz  = data_3hz.s;  % 7 × 110
        partial_corr_3hz = data_3hz.p;  % 7 × 110
        
        s_liking_3hz      = [s_liking_3hz;      simple_corr_3hz(1, :)];
        s_valence_3hz     = [s_valence_3hz;     simple_corr_3hz(2, :)];
        s_arousal_3hz     = [s_arousal_3hz;     simple_corr_3hz(3, :)];
        s_complexity_3hz  = [s_complexity_3hz;  simple_corr_3hz(4, :)];
        s_familiarity_3hz = [s_familiarity_3hz; simple_corr_3hz(5, :)];
        s_style_3hz       = [s_style_3hz;       simple_corr_3hz(6, :)];
        s_category_3hz    = [s_category_3hz;    simple_corr_3hz(7, :)];
        
        p_liking_3hz      = [p_liking_3hz;      partial_corr_3hz(1, :)];
        p_valence_3hz     = [p_valence_3hz;     partial_corr_3hz(2, :)];
        p_arousal_3hz     = [p_arousal_3hz;     partial_corr_3hz(3, :)];
        p_complexity_3hz  = [p_complexity_3hz;  partial_corr_3hz(4, :)];
        p_familiarity_3hz = [p_familiarity_3hz; partial_corr_3hz(5, :)];
        p_style_3hz       = [p_style_3hz;       partial_corr_3hz(6, :)];
        p_category_3hz    = [p_category_3hz;    partial_corr_3hz(7, :)];
    else
        fprintf('File not found: %s\n', filename_3hz);
    end
end

%%------------------------ Load Time Vector -------------------------------
% Load the time vector (1 × 110) for plotting
time_data = load('D:\Dropbox\Internship\tasks_may2024\time.mat');
time_vec  = time_data.time;

%%------------------------- Plot Settings ---------------------------------
% Feature names and corresponding colors for plotting
variables = {'Liking', 'Valence', 'Arousal', 'Complexity', 'Familiarity', 'Style', 'Category'};
feature_colors = containers.Map(variables, { ...
    [0.8500, 0.3250, 0.0980], ...   % Liking
    [0,      0.4470, 0.7410], ...   % Valence
    [0.4660, 0.6740, 0.1880], ...   % Arousal
    [0.9290, 0.6940, 0.1250], ...   % Complexity
    [0.4940, 0.1840, 0.5560], ...   % Familiarity
    [0.3010, 0.7450, 0.9330], ...   % Style
    [0.6350, 0.0780, 0.1840]  ...   % Category
});
lineWidths = [3, 2, 1];  % [highlighted-feature width, other-feature width, marker width]
highlightFeature = 'Liking';  % Feature to highlight in plots

%%==================== Plot: Partial Correlation ==========================
figure;
tiledlayout(1, 2, 'TileSpacing', 'Compact', 'Padding', 'Compact');

% ------- Left subplot: 20 Hz Partial Correlations -------
nexttile;
hold on;
legend_handles = [];

% Plot average_partial for each feature over time
for idx = 1:length(variables)
    feature_name = variables{idx};
    color_full   = feature_colors(feature_name);
    color_faded  = color_full * 0.5 + 0.5;
    
    if strcmp(feature_name, highlightFeature)
        h = plot(time_vec, average_partial(idx, :), ...
            'Color', color_full, ...
            'LineStyle', '-', ...
            'LineWidth', lineWidths(2), ...
            'DisplayName', feature_name);
    else
        h = plot(time_vec, average_partial(idx, :), ...
            'Color', color_faded, ...
            'LineStyle', '-', ...
            'LineWidth', lineWidths(3), ...
            'DisplayName', feature_name);
    end
    
    legend_handles = [legend_handles, h];
end

% Shade stimulus-on period (0 to 50 ms)
p = fill([0, 0.050, 0.050, 0], [-0.0105, -0.0105, -0.0095, -0.0095], ...
    [0.5, 0.5, 0.5], 'EdgeColor', 'k');
p.LineWidth = 2;

% Dummy patch object for legend ('Stimulus on')
h_stim = patch(NaN, NaN, [0.5, 0.5, 0.5], 'DisplayName', 'Stimulus on', 'LineWidth', 1);
legend_handles = [legend_handles, h_stim];

% Perform one-sample, right-tailed t-tests per feature across participants
feature_data_names = { 'p_liking',   'p_valence',   'p_arousal',   ...
                       'p_complexity','p_familiarity','p_style', ...
                       'p_category' };
p_adj_matrix = [];

for i = 1:length(feature_data_names)
    data_matrix = eval(feature_data_names{i});  % e.g., p_liking is Nparticipants × 110
    [~, p_vals]  = ttest(data_matrix, 0, 'tail', 'right');
    [~, ~, ~, adj_p] = fdr_bh(p_vals);  % Benjamini-Hochberg correction
    
    p_adj_matrix = [p_adj_matrix; adj_p];
    
    % Determine significant timepoints (p < 0.05 and time > 0)
    sig_timepoints = (adj_p < 0.05) & (time_vec > 0);
    offset = (i - 1) * 0.002; 
    ylims  = [min(average_partial(:)) + offset, max(average_partial(:)) + offset];
    
    if any(sig_timepoints)
        if strcmp(feature_data_names{i}, 'p_liking')
            marker_color = feature_colors(variables{i});
        else
            marker_color = feature_colors(variables{i}) * 0.5 + 0.5;
        end
        
        hsig = plot(time_vec(sig_timepoints), ylims(2) + 0.002, 's', ...
            'MarkerSize', 4, ...
            'MarkerEdgeColor', marker_color, ...
            'MarkerFaceColor', marker_color, ...
            'LineWidth', 1);
    end
end

% Customize axes
xline(0, '--');
yline(0, '--');
xlim([min(time_vec), max(time_vec)]);
xlabel('Time (s)');
ylabel('Partial Correlation');
title('50 ms Partial Correlation (20 Hz)');
legend(legend_handles, 'Location', 'best');
grid off;
hold off;

% ------- Right subplot: 3 Hz Partial Correlations -------
nexttile;
hold on;
legend_handles_3hz = [];

for idx = 1:length(variables)
    feature_name = variables{idx};
    color_full   = feature_colors(feature_name);
    color_faded  = color_full * 0.5 + 0.5;
    
    if strcmp(feature_name, highlightFeature)
        h = plot(time_vec, average_partial_3hz(idx, :), ...
            'Color', color_full, ...
            'LineStyle', '-', ...
            'LineWidth', lineWidths(2), ...
            'DisplayName', feature_name);
    else
        h = plot(time_vec, average_partial_3hz(idx, :), ...
            'Color', color_faded, ...
            'LineStyle', '-', ...
            'LineWidth', lineWidths(3), ...
            'DisplayName', feature_name);
    end
    
    legend_handles_3hz = [legend_handles_3hz, h];
end

% Shade stimulus-on (0 to 150 ms)
p = fill([0, 0.150, 0.150, 0], [-0.0105, -0.0105, -0.0095, -0.0095], ...
    [0.5, 0.5, 0.5], 'EdgeColor', 'k');
p.LineWidth = 2;

h_stim = patch(NaN, NaN, [0.5, 0.5, 0.5], 'DisplayName', 'Stimulus on', 'LineWidth', 1);
legend_handles_3hz = [legend_handles_3hz, h_stim];

feature_data_names_3hz = { 'p_liking_3hz',   'p_valence_3hz',   'p_arousal_3hz',   ...
                           'p_complexity_3hz','p_familiarity_3hz','p_style_3hz', ...
                           'p_category_3hz' };
p_adj_matrix_3hz = [];

for i = 1:length(feature_data_names_3hz)
    data_matrix_3hz = eval(feature_data_names_3hz{i});
    [~, p_vals_3hz]    = ttest(data_matrix_3hz, 0, 'tail', 'right');
    [~, ~, ~, adj_p_3hz] = fdr_bh(p_vals_3hz);
    
    p_adj_matrix_3hz = [p_adj_matrix_3hz; adj_p_3hz];
    
    sig_timepoints_3hz = (adj_p_3hz < 0.05) & (time_vec > 0);
    
    if any(sig_timepoints_3hz)
        if strcmp(feature_data_names_3hz{i}, 'p_liking_3hz')
            marker_color = feature_colors(variables{i});
        else
            marker_color = feature_colors(variables{i}) * 0.5 + 0.5;
        end
        
        offset = (i - 1) * 0.002;
        ylims  = [min(average_partial_3hz(:)) + offset, max(average_partial_3hz(:)) + offset];
        
        hsig = plot(time_vec(sig_timepoints_3hz), ylims(2) + 0.002, 's', ...
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
ylabel('Partial Correlation');
title('150 ms Partial Correlation (3 Hz)');
legend(legend_handles_3hz, 'Location', 'best');
grid off;
hold off;

% Save the partial correlation figure
% plotpath_partial = fullfile('D:\Dropbox\Internship\tasks_may2024\Results\average\', ...
%                             'Partial_correlation.png');
% set(gcf, 'Position', [100, 100, 1200, 600]);
% exportgraphics(gcf, plotpath_partial, 'Resolution', 300);

%%===================== Plot: Simple Correlation ==========================
figure;
tiledlayout(1, 2, 'TileSpacing', 'Compact', 'Padding', 'Compact');

% ------- Left subplot: 20 Hz Simple Correlations -------
nexttile;
hold on;
legend_handles_simple = [];

for idx = 1:length(variables)
    feature_name = variables{idx};
    color_full   = feature_colors(feature_name);
    color_faded  = color_full * 0.5 + 0.5;
    
    if strcmp(feature_name, highlightFeature)
        h = plot(time_vec, average_simple(idx, :), ...
            'Color', color_full, ...
            'LineStyle', '-', ...
            'LineWidth', lineWidths(2), ...
            'DisplayName', feature_name);
    else
        h = plot(time_vec, average_simple(idx, :), ...
            'Color', color_faded, ...
            'LineStyle', '-', ...
            'LineWidth', lineWidths(3), ...
            'DisplayName', feature_name);
    end
    
    legend_handles_simple = [legend_handles_simple, h];
end

% Shade stimulus-on (0 to 50 ms) slightly below the min value
min_simple = min(average_simple(:));
p = fill([0, 0.050, 0.050, 0], [min_simple-0.003, min_simple-0.003, -0.005, -0.005], ...
    [0.5, 0.5, 0.5], 'EdgeColor', 'k');
p.LineWidth = 2;

h_stim = patch(NaN, NaN, [0.5, 0.5, 0.5], 'DisplayName', 'Stimulus on', 'LineWidth', 1);
legend_handles_simple = [legend_handles_simple, h_stim];

feature_data_names_simple = { 's_liking',   's_valence',   's_arousal',   ...
                              's_complexity','s_familiarity','s_style', ...
                              's_category' };
s_adj_matrix = [];

for i = 1:length(feature_data_names_simple)
    data_matrix_simple = eval(feature_data_names_simple{i});
    [~, p_vals_simple]    = ttest(data_matrix_simple, 0, 'tail', 'right');
    [~, ~, ~, adj_p_simple] = fdr_bh(p_vals_simple);
    s_adj_matrix = [s_adj_matrix; adj_p_simple];
    
    sig_timepoints_simple = (adj_p_simple < 0.05) & (time_vec > 0);
    if any(sig_timepoints_simple)
        if strcmp(feature_data_names_simple{i}, 's_liking')
            marker_color = feature_colors(variables{i});
        else
            marker_color = feature_colors(variables{i}) * 0.5 + 0.5;
        end
        
        offset = (i - 1) * 0.002;
        ylims  = [min_simple + offset, max(average_simple(:)) + offset];
        
        hsig = plot(time_vec(sig_timepoints_simple), ylims(2) + 0.002, 's', ...
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
ylabel('Simple Correlation');
title('50 ms Simple Correlation (20 Hz)');
legend(legend_handles_simple, 'Location', 'best');
grid off;
hold off;

% ------- Right subplot: 3 Hz Simple Correlations -------
nexttile;
hold on;
legend_handles_simple_3hz = [];

for idx = 1:length(variables)
    feature_name = variables{idx};
    color_full   = feature_colors(feature_name);
    color_faded  = color_full * 0.5 + 0.5;
    
    if strcmp(feature_name, highlightFeature)
        h = plot(time_vec, average_simple_3hz(idx, :), ...
            'Color', color_full, ...
            'LineStyle', '-', ...
            'LineWidth', lineWidths(2), ...
            'DisplayName', feature_name);
    else
        h = plot(time_vec, average_simple_3hz(idx, :), ...
            'Color', color_faded, ...
            'LineStyle', '-', ...
            'LineWidth', lineWidths(3), ...
            'DisplayName', feature_name);
    end
    
    legend_handles_simple_3hz = [legend_handles_simple_3hz, h];
end

min_simple_3hz = min(average_simple_3hz(:));
p = fill([0, 0.150, 0.150, 0], [min_simple_3hz-0.003, min_simple_3hz-0.003, -0.005, -0.005], ...
    [0.5, 0.5, 0.5], 'EdgeColor', 'k');
p.LineWidth = 2;

h_stim = patch(NaN, NaN, [0.5, 0.5, 0.5], 'DisplayName', 'Stimulus on', 'LineWidth', 1);
legend_handles_simple_3hz = [legend_handles_simple_3hz, h_stim];

s_adj_matrix_3hz = [];

for i = 1:length(feature_data_names_simple)
    data_matrix_simple_3hz = eval([feature_data_names_simple{i}, '_3hz']);
    [~, p_vals_simple_3hz]     = ttest(data_matrix_simple_3hz, 0, 'tail', 'right');
    [~, ~, ~, adj_p_simple_3hz] = fdr_bh(p_vals_simple_3hz);
    s_adj_matrix_3hz = [s_adj_matrix_3hz; adj_p_simple_3hz];
    
    sig_timepoints_simple_3hz = (adj_p_simple_3hz < 0.05) & (time_vec > 0);
    
    if any(sig_timepoints_simple_3hz)
        if strcmp(feature_data_names_simple{i}, 's_liking')
            marker_color = feature_colors(variables{i});
        else
            marker_color = feature_colors(variables{i}) * 0.5 + 0.5;
        end
        
        offset = (i - 1) * 0.002;
        ylims  = [min_simple_3hz + offset, max(average_simple_3hz(:)) + offset];
        
        hsig = plot(time_vec(sig_timepoints_simple_3hz), ylims(2) + 0.002, 's', ...
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
ylabel('Simple Correlation');
title('150 ms Simple Correlation (3 Hz)');
legend(legend_handles_simple_3hz, 'Location', 'best');
grid off;
hold off;

% Save the simple correlation figure
% plotpath_simple = fullfile('D:\Dropbox\Internship\tasks_may2024\Results\average\', ...
%                            'Simple_Correlation.png');
% set(gcf, 'Position', [100, 100, 1200, 600]);
% exportgraphics(gcf, plotpath_simple, 'Resolution', 300);

%%============ Plot: Before vs After (20 Hz, Partial & Simple) =============
% Partial correlations before vs. after (20 Hz)
figure;
tiledlayout(1, 2, 'TileSpacing', 'Compact', 'Padding', 'Compact');

% --- Partial Correlation (Before) ---
nexttile;
hold on;
legend_handles_bef = [];

for idx = 1:length(variables)
    feature_name = variables{idx};
    color_full   = feature_colors(feature_name);
    
    h = plot(time_vec, average_partial_before(idx, :), ...
        'Color', color_full, ...
        'LineStyle', '-', ...
        'LineWidth', lineWidths(2), ...
        'DisplayName', feature_name);
    
    legend_handles_bef = [legend_handles_bef, h];
end

% Shade stimulus-on (0 to 50 ms)
p = fill([0, 0.050, 0.050, 0], [-0.0035, -0.0035, -0.005, -0.005], ...
    [0.5, 0.5, 0.5], 'EdgeColor', 'k');
p.LineWidth = 2;

h_stim = patch(NaN, NaN, [0.5, 0.5, 0.5], 'DisplayName', 'Stimulus on', 'LineWidth', 1);
legend_handles_bef = [legend_handles_bef, h_stim];

xline(0, '--');
yline(0, '--');
xlim([min(time_vec), max(time_vec)]);
xlabel('Time (s)');
ylabel('Partial Correlation');
title('50 ms Partial Corr. Before (20 Hz)');
legend(legend_handles_bef, 'Location', 'best');
grid off;
hold off;

% --- Partial Correlation (After) ---
nexttile;
hold on;
legend_handles_aft = [];

for idx = 1:length(variables)
    feature_name = variables{idx};
    color_full   = feature_colors(feature_name);
    
    h = plot(time_vec, average_partial_after(idx, :), ...
        'Color', color_full, ...
        'LineStyle', '-', ...
        'LineWidth', lineWidths(2), ...
        'DisplayName', feature_name);
    
    legend_handles_aft = [legend_handles_aft, h];
end

% Shade stimulus-on (0 to 150 ms)
p = fill([0, 0.150, 0.150, 0], [-0.0035, -0.0035, -0.005, -0.005], ...
    [0.5, 0.5, 0.5], 'EdgeColor', 'k');
p.LineWidth = 2;

h_stim = patch(NaN, NaN, [0.5, 0.5, 0.5], 'DisplayName', 'Stimulus on', 'LineWidth', 1);
legend_handles_aft = [legend_handles_aft, h_stim];

xline(0, '--');
yline(0, '--');
xlim([min(time_vec), max(time_vec)]);
xlabel('Time (s)');
ylabel('Partial Correlation');
title('50 ms Partial Corr. After (20 Hz)');
legend(legend_handles_aft, 'Location', 'best');
grid off;
hold off;

% Save before/after partial correlations figure
% plotpath_befaft_partial = fullfile('D:\Dropbox\Internship\tasks_may2024\Results\average\', ...
%                                    'Partial_correlation_before_after.png');
% set(gcf, 'Position', [100, 100, 1200, 600]);
% exportgraphics(gcf, plotpath_befaft_partial, 'Resolution', 300);

%%-------------------- Simple Correlation Before/After --------------------
figure;
tiledlayout(1, 2, 'TileSpacing', 'Compact', 'Padding', 'Compact');

% --- Simple Correlation (Before) ---
nexttile;
hold on;
legend_handles_bef = [];

for idx = 1:length(variables)
    feature_name = variables{idx};
    color_full   = feature_colors(feature_name);
    
    h = plot(time_vec, average_simple_before(idx, :), ...
        'Color', color_full, ...
        'LineStyle', '-', ...
        'LineWidth', lineWidths(2), ...
        'DisplayName', feature_name);
    
    legend_handles_bef = [legend_handles_bef, h];
end

% Shade stimulus-on (0 to 50 ms)
p = fill([0, 0.050, 0.050, 0], [-0.0035, -0.0035, -0.005, -0.005], ...
    [0.5, 0.5, 0.5], 'EdgeColor', 'k');
p.LineWidth = 2;

h_stim = patch(NaN, NaN, [0.5, 0.5, 0.5], 'DisplayName', 'Stimulus on', 'LineWidth', 1);
legend_handles_bef = [legend_handles_bef, h_stim];

xline(0, '--');
yline(0, '--');
xlim([min(time_vec), max(time_vec)]);
xlabel('Time (s)');
ylabel('Simple Correlation');
title('50 ms Simple Corr. Before (20 Hz)');
legend(legend_handles_bef, 'Location', 'best');
grid off;
hold off;

% --- Simple Correlation (After) ---
nexttile;
hold on;
legend_handles_aft = [];

for idx = 1:length(variables)
    feature_name = variables{idx};
    color_full   = feature_colors(feature_name);
    
    h = plot(time_vec, average_simple_after(idx, :), ...
        'Color', color_full, ...
        'LineStyle', '-', ...
        'LineWidth', lineWidths(2), ...
        'DisplayName', feature_name);
    
    legend_handles_aft = [legend_handles_aft, h];
end

% Shade stimulus-on (0 to 150 ms)
p = fill([0, 0.150, 0.150, 0], [-0.0035, -0.0035, -0.005, -0.005], ...
    [0.5, 0.5, 0.5], 'EdgeColor', 'k');
p.LineWidth = 2;

h_stim = patch(NaN, NaN, [0.5, 0.5, 0.5], 'DisplayName', 'Stimulus on', 'LineWidth', 1);
legend_handles_aft = [legend_handles_aft, h_stim];

xline(0, '--');
yline(0, '--');
xlim([min(time_vec), max(time_vec)]);
xlabel('Time (s)');
ylabel('Simple Correlation');
title('50 ms Simple Corr. After (20 Hz)');
legend(legend_handles_aft, 'Location', 'best');
grid off;
hold off;

% Save before/after simple correlations figure
% plotpath_befaft_simple = fullfile('D:\Dropbox\Internship\tasks_may2024\Results\average\', ...
%                                   'Simple_correlation_before_after.png');
% set(gcf, 'Position', [100, 100, 1200, 600]);
% exportgraphics(gcf, plotpath_befaft_simple, 'Resolution', 300);
