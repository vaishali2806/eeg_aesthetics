%--------------------------------------------------------------------------
% This script:
%  1. Loads VAPS metadata to build binary style and category model RDMs.
%  2. Averages precomputed EEG RDMs across participants (1–39) at 3 Hz.
%  3. Computes Spearman‐based RSA (via corr + Fisher z) between the averaged
%     EEG RDM and each style/category model RDM at each time bin.
%  4. Plots the timecourses of RSA correlations for styles and categories.
%  5. Exports figures to the Results/average folder.
%--------------------------------------------------------------------------

clear all;
clc;

%%======================= LOAD VAPS METADATA ===============================
% Load database_VAPS.mat which contains dataVaps.data (999×N_features)
dataVaps = load("D:\Dropbox\Internship\tasks_may2024\database_VAPS.mat");
style_col = dataVaps.data(:, 3);
cat_col = dataVaps.data(:, 4);

% GROUP STYLES IN 3 CATEGORIES
group1 = [1, 2, 3, 4, 9];
group2 = [5, 6, 7, 10];
group3 = [8, 11, 12, 13];
new_style = zeros(size(style_col));

new_style(ismember(style_col, group1)) = 1;
new_style(ismember(style_col, group2)) = 2;
new_style(ismember(style_col, group3)) = 3;

dataVaps.data(:, 3) = new_style;
%%
%%===================== BUILD BINARY MODEL RDMs ===========================

style_length = length(new_style);
cat_length = length(cat_col);
cat1_matrix = zeros(cat_length);
cat2_matrix = zeros(cat_length);
cat3_matrix = zeros(cat_length);
cat4_matrix = zeros(cat_length);
cat5_matrix = zeros(cat_length);
style1_matrix = zeros(style_length);
style2_matrix = zeros(style_length);
style3_matrix = zeros(style_length);


is_cat1 = (cat_col == 1);
is_cat2 = (cat_col == 2);
is_cat3 = (cat_col == 3);
is_cat4 = (cat_col == 4);
is_cat5 = (cat_col == 5);
is_style1 = (new_style == 1);
is_style2 = (new_style == 2);
is_style3 = (new_style == 3);


cat1_matrix = 1 - double(is_cat1 * is_cat1');
cat2_matrix = 1 - double(is_cat2 * is_cat2');
cat3_matrix = 1 - double(is_cat3 * is_cat3');
cat4_matrix = 1 - double(is_cat4 * is_cat4');
cat5_matrix = 1 - double(is_cat5 * is_cat5');
style1_matrix = 1 - double(is_style1 * is_style1');
style2_matrix = 1 - double(is_style2 * is_style2');
style3_matrix = 1 - double(is_style3 * is_style3');

%%
%%====================== AVERAGE EEG RDMs ACROSS SUBJECTS FOR 3 HZ ==================
n_participants = 39;
n_timebins = 110;
matrix_size = 999;

eeg_avg = zeros(matrix_size, matrix_size, n_timebins);  
valid_counts = 0;

for p = 1:n_participants
    filename = fullfile(strcat("G:\.shortcut-targets-by-id\1uGZWB6XgxS7bqYQVdYfB_UgkXp8VePvc\Thesis_VG\Participants analysis\eeg_aesthetics_", string(p)), ...
        "3Hz", strcat("ds_corr_", string(p)));
    if exist("filename", "file")
        loaded = load(filename);
        
        eeg_avg = eeg_avg + loaded.ds_corr;  
        valid_counts = valid_counts + 1;
        
        disp(['Loaded and added participant ' num2str(p)])
    else
        disp(['participant ' num2str(p) ' not found'])
    end
end

% Compute the mean EEG RDM across participants
eeg_avg = eeg_avg / valid_counts;  % result: [999 x 999 x 110]
%%
model1_vec = style1_matrix(:);
model2_vec = style2_matrix(:);
model3_vec = style3_matrix(:);
model1_cat = cat1_matrix(:);
model2_cat = cat2_matrix(:);
model3_cat = cat3_matrix(:);
model4_cat = cat4_matrix(:);
model5_cat = cat5_matrix(:);

mask = ~eye(matrix_size);
valid_idx = mask(:);

corrs_style1 = zeros(n_timebins,1);
corrs_style2 = zeros(n_timebins,1);
corrs_style3 = zeros(n_timebins,1);
corrs_cat1 = zeros(n_timebins, 1);
corrs_cat2 = zeros(n_timebins, 1);
corrs_cat3 = zeros(n_timebins, 1);
corrs_cat4 = zeros(n_timebins, 1);
corrs_cat5 = zeros(n_timebins, 1);

for t = 1:n_timebins
    eeg_vec = eeg_avg(:,:,t);
    eeg_vec = eeg_vec(:);

    corrs_style1(t) = corr(eeg_vec(valid_idx), model1_vec(valid_idx));
    corrs_style2(t) = corr(eeg_vec(valid_idx), model2_vec(valid_idx));
    corrs_style3(t) = corr(eeg_vec(valid_idx), model3_vec(valid_idx));
    corrs_cat1(t) = corr(eeg_vec(valid_idx), model1_cat(valid_idx));
    corrs_cat2(t) = corr(eeg_vec(valid_idx), model2_cat(valid_idx));
    corrs_cat3(t) = corr(eeg_vec(valid_idx), model3_cat(valid_idx));
    corrs_cat4(t) = corr(eeg_vec(valid_idx), model4_cat(valid_idx));
    corrs_cat5(t) = corr(eeg_vec(valid_idx), model5_cat(valid_idx));

     if (rem(t,10) ==0)
            disp(['time bin ' num2str(t) ' processed'])
    end 
end
%%
%%======================== PLOT: STYLE‐MODEL RSA ===========================
time = 1:n_timebins;
num_style = 3;
num_cat = 5;
colors = lines(5);
legend_styles = {'Classical Tendencies', 'Romantic/Emotive Tendencies', 'Modernist Tendencies'};
figure;
hold on;
for i = 1: num_style
    style_data = eval(strcat("corrs_style", string(i)));
    plot(time, style_data, 'Color', colors(i, :), 'LineWidth', 1);

end

xline(0, '--'); 
yline(0, '--');
xlabel('Time Bin');
ylabel('RSA Correlation');
legend(legend_styles);
title('150ms RSA: Averaged EEG RDMs vs Style Models');
grid off;
hold off;
plotpath = fullfile("D:\Dropbox\Internship\tasks_may2024\Results\average", 'corr_style_liking_modeling.png');
set(gcf, 'Position', [100, 100, 1200, 600]); 
exportgraphics(gcf, plotpath, 'Resolution', 300); 
%%
%%======================== PLOT: CATEGORY‐MODEL RSA ===========================
legend_cat = {"Scenes", "Portrait", "Landscape", "Still life", "Toward Abstraction"};
figure;
hold on;
for i = 1: num_cat
    cat_data = eval(strcat("corrs_cat", string(i)));
    plot(time, cat_data, 'Color', colors(i, :), 'LineWidth', 1);

end

xline(0, '--'); 
yline(0, '--');
xlabel('Time Bin');
ylabel('RSA Correlation');
legend(legend_cat);
title('150ms RSA: Averaged EEG RDMs vs Category Models');
grid off;
hold off;

plotpath = fullfile("D:\Dropbox\Internship\tasks_may2024\Results\average", 'corr_category_liking_modeling.png');
set(gcf, 'Position', [100, 100, 1200, 600]); 
exportgraphics(gcf, plotpath, 'Resolution', 300); 