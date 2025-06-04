%%
%--------------------------------------------------------------------------
% This script performs the following steps for a single subject:
%  1. Prompts for subject number and creates an output folder for 3 Hz results.
%  2. Loads behavioral log (RSVP triggers) and maps each trial to its VAPS ratings.
%  3. Uses FieldTrip to:
%       • Define trials around stimulus events,
%       • Preprocess EEG (bandstop, re-reference, downsample),
%       • Identify and remove flat‐channels and artifactual ICA components,
%       • Epoch, baseline‐correct, and save the cleaned timelocked data.
%  4. Builds a CoSMoMVPA dataset, runs PCA (99 % variance) per time‐bin, averages samples,
%     and computes a Spearman distance matrix (ds_corr) for each time‐bin.
%  5. Computes “simple” (Spearman‐based) and “partial” correlations between:
%       • Each VAPS dimension and the EEG RDM,
%       • Binary style and category predictors and the EEG RDM.
%  6. Plots and saves timecourses of all seven “simple” correlations.
%  7. Plots and saves timecourses of all seven “partial” correlations.
%  8. Separately computes and plots “style vs. liking” and “category vs. liking”
%     correlations for each style/category label.
%--------------------------------------------------------------------------

clc
close all
clear all
%%========================== ADD PATHS ====================================
% Add FieldTrip and CoSMoMVPA to MATLAB path
addpath(genpath('C:\Masters Cyber\Master Thesis\toolboxes\CoSMoMVPA-master'));
ft_defaults();  % initialize FieldTrip defaults

%%====================== INPUT AND OUTPUT FOLDER ==========================
%=========================== USER SETUP =============================
freqChoice = input('Which frequency do you want to process? Enter 20 or 3: ');
if freqChoice == 20
    freqStr = '20Hz';
elseif freqChoice == 3
    freqStr = '3Hz';
else
    error('Invalid input. Please enter 20 or 3.');
end
% Prompt for subject number
participant_num = input('Enter Subject Number: ');
% 4) Build the base folder for this subject
baseFolder = sprintf('C:\\Masters Cyber\\Master Thesis\\Participants analysis\\eeg_aesthetics_%d', participant_num);

% 5) Inside that, create (or check) a subfolder named either "20Hz" or "3Hz"
foldername = fullfile(baseFolder, freqStr);
if ~exist(foldername, 'dir')
    mkdir(foldername);
    fprintf('Created folder: %s\n', foldername);
else
    fprintf('Using existing folder: %s\n', foldername);
end
%%========================== LOAD LOG FILE =================================
% Load the RSVP log file for this subject, to map trial indices to VAPS ratings
if strcmp(freqStr, '20Hz')
    logfile = load(strcat('C:\Masters Cyber\Master Thesis\Data-VG\log\RSVP_eeg_s',string(participant_num), 'f.mat'));
else
    logfile = load(strcat('C:\Masters Cyber\Master Thesis\Data-VG\log\RSVP_eeg_s',string(participant_num), '.mat'));
end
%%
% Load entire VAPS database (subjects × features)
data_VAPS = load('C:\Masters Cyber\Master Thesis\database_VAPS.mat');
 %% 
 %%======================== ASSIGN VAPS TO TRIALS ===========================
% Preallocate a matrix to hold each trial’s 1×N‐feature VAPS vector
% (we will fill only when trigger == 1)

for ii = 1: length(logfile.dat.new_stim)
    % Only proceed if this trial’s trigger is “1” (stimulus onset)
    if logfile.dat.trigs(ii) ==1
        image_presented = logfile.dat.new_stim{ii};
        
        % Extract the numeric ID from the image filename
        pattern = '\d+';
        number = regexp(image_presented, pattern, 'match');
        
        % Convert the extracted number from cell to double
        if ~isempty(number)
            extracted_number = str2double(number{1});
        else
            extracted_number = [];
        end
        
        % Map that image’s row in data_VAPS.data into trials_data      
        trials_data(ii,:) = data_VAPS.data(extracted_number,:);
    end
end
% 

%%
%%===================== DEFINE FIELDTRIP TRIALS ===========================
% Specify the raw .eeg file for this subject
if strcmp(freqStr, '20Hz')
    fileName = (strcat('C:\Masters Cyber\Master Thesis\Data-VG\eeg\eeg_aesthetics_00', string(participant_num),'f.eeg'));
else
    if (participant_num < 10)
        fileName = (strcat('C:\Masters Cyber\Master Thesis\Data-VG\eeg\eeg_aesthetics_000', string(participant_num),'.eeg'));
    else
        fileName = (strcat('C:\Masters Cyber\Master Thesis\Data-VG\eeg\eeg_aesthetics_00', string(participant_num),'.eeg'));
    end
end
%Define Events
cfg=[]; 
cfg.dataset=char(fileName);
cfg.trialdef.eventtype= 'Stimulus';

cfg.trialdef.prestim=0.1;
cfg.trialdef.poststim=1;
cfg=ft_definetrial(cfg);
%%
%%============================ PREPROCESSING ================================
% Only keep EEG channels; apply bandstop at 50 Hz, de-mean, re-reference to FCz
cfg.channel='eeg';
cfg.hpfilter='no';
cfg.lpfilter='no';
cfg.bsfilter='yes';
cfg.bsfreq=[48,52];
cfg.demean='yes';
%cfg.baselinewindow=[-0.25,0];
cfg.reref='yes';
cfg.refchannel='FCz';
data=ft_preprocessing(cfg);
%%
if(size(data.trialinfo, 1) - size(trials_data,1) > 0)
    
    valid_indices = ismember(data.trialinfo, [1, 2]);
    data.trialinfo = data.trialinfo(valid_indices, :);
    data.time = data.time(:, valid_indices);
    data.trial = data.trial(:, valid_indices);
    data.sampleinfo = data.sampleinfo(valid_indices, :);
end
data.trialinfo  = [data.trialinfo trials_data];
%%
%%========================= DOWNSAMPLE TO 100 Hz ============================
cfg.resamplefs = 100; 
data_selected = ft_resampledata(cfg,data);
%%
%%======================= REMOVE FLAT (ZERO‐VAR) CHANNELS ==================
zero_var_channels = [];
for i = 1:length(data_selected.label)
    channel_variances = zeros(1, length(data_selected.trial));
    for j = 1:length(data_selected.trial)
        channel_data = data_selected.trial{j}(i, :);
        channel_variances(j) = var(channel_data);
    end
    % If variance is zero across all trials, mark channel for removal
    if all(channel_variances == 0)
        zero_var_channels = [zero_var_channels; data_selected.label(i)];
    end
end

% Display the channels with zero variance (flat signals)
disp('Channels with zero frequency (zero variance):');
disp(zero_var_channels);
%%
% Remove channels with zero variance from the data
cfg = [];
cfg.channel = setdiff(data_selected.label, zero_var_channels); 
clean_data = ft_selectdata(cfg, data_selected);
%%
%%========================= VISUAL ARTIFACT REJECTION =====================
cfg=[];
cfg.showlabel='yes';
cfg.method='summary'; 
cfg.keepchannel='no';
data_clean=ft_rejectvisual(cfg,clean_data);
%%
%%========================== RUN ICA AND REMOVE ============================
cfg = [];
cfg.method = 'runica';
cfg.numcomponent = 30; 
comp = ft_componentanalysis(cfg, data_clean);
%%
% Plot topographies for the first 30 components
cfg.component = 1:30; 
cfg.layout = 'easycap-M1.txt'; 
cfg.marker = 'off';
ft_topoplotIC(cfg, comp);
%%
% Inspect component timecourses as needed:
cfg = [];
cfg.viewmode = 'component';
cfg.channel = 4;
cfg.layout    = 'easycap-M1.txt';
ft_databrowser(cfg, comp);
%%
% Manually decide which component(s) to remove
cfg = [];
cfg.component = 3;
data_clean = ft_rejectcomponent(cfg, comp, data_clean);
%%
%%==================== EPOCH SELECTION (CONDITION == 1) ====================
% Keep only trials where trialinfo(:,1) == 1 (condition code for “stimulus present”)
cfg= [];
cfg.trials = find(data_clean.trialinfo(:,1) ==1);
data_subset = ft_selectdata(cfg,data_clean);


% Compute timelock (ERP) per trial
cfg=[];
cfg.keeptrials = 'yes';
data_subset = ft_timelockanalysis(cfg,data_subset);
%% 
%%========================== BASELINE CORRECTION ===========================
cfg=[];
cfg.baseline = [-0.05 0];
data_subset_bs = ft_timelockbaseline(cfg,data_subset);

% data_subset_bs_filename = fullfile(foldername, strcat('data_subset_bs_', string(participant_num)));
% save(data_subset_bs_filename, 'data_subset_bs');
%%
%%===================== BUILD CoSMoMVPA DATASET ============================
% Convert the FieldTrip timelocked data into a CoSMoMVPA dataset
ds=cosmo_meeg_dataset(data_subset_bs);
ds.sa.targets = data_subset.trialinfo(:,2);
ds.sa.chunks = ones(length(data_subset.trialinfo),1);
%% 
%%================== PCA + AVERAGE PER TIME BIN + RDM =====================
for n_bin = 1:size(data_subset_bs.trial,3)
    % Slice dataset to this time bin (feature dimension = channels)
    ds_t = cosmo_slice(ds,ismember(ds.fa.time,n_bin),2);
    
    % Run PCA, keep 99 % variance, max 100 000 features
    ds_pca = cosmo_map_pca(ds_t, 'pca_explained_ratio', 0.99, 'max_feature_count', 100000);
    ds_pca.sa.targets = data_subset.trialinfo(:,2);
    ds_pca.sa.chunks = ones(length(data_subset.trialinfo),1);

    % Average all samples (across trials) into one “average sample” per target
    ds_avg_samples = cosmo_average_samples(ds_pca);
    
    % Compute Spearman distance: 1 − corr(matrix, 'Spearman')
    ds_corr(:,:,n_bin) = 1 - cosmo_corr(ds_avg_samples.samples','Spearman');
    if (rem(n_bin,50) ==0)
        disp(['time bin ' num2str(n_bin) ' processed'])
    end 
end

% ds_corr_filename = fullfile(foldername, strcat('ds_corr_', string(participant_num)));
% save(ds_corr_filename, 'ds_corr');
%% 
%%==================== BUILD PREDICTOR MATRICES ============================
% “pred_behav” will store absolute differences between every pair of stimuli
% for each VAPS dimension of interest: Liking (col 8), Valence (10), Arousal (12),
% Complexity (14), Familiarity (16).
cond_colms = [8, 10, 12, 14,16];
for n_cond = 1:5
for ii = 1:999 
    pred_val_1 = data_VAPS.data(ii,cond_colms(n_cond));
    for jj = 1:999
        pred_val_2 = data_VAPS.data(jj,cond_colms(n_cond));
        diff_pred(ii,jj) = abs(pred_val_1 - pred_val_2);
    end 
end
% Convert full 999×999 to condensed (lower triangle) as a column vector
pred_behav(:,n_cond) = squareform(diff_pred)';
end 
%% 
%%================== SIMPLE (SPEARMAN) CORRELATION ========================
numRowsPredBehav = size(pred_behav, 1);
numColsPredBehav = size(pred_behav, 2);
numTimeBins = size(ds_corr, 3);

corr_behav_eeg = zeros(numColsPredBehav, numTimeBins);

for ii =1: numTimeBins
    % Extract the EEG distance matrix for time bin t & vectorize
    eeg_pred = squareform(squeeze(ds_corr(:,:,ii)))';
    temp_corr = zeros(numColsPredBehav, 1);
    for kk = 1:numColsPredBehav
        % Compute Spearman correlation, then Fisher‐z (atanh) for normality
        temp_corr(kk) = atanh(cosmo_corr(pred_behav(:, kk), eeg_pred, 'Spearman'));
    end
    corr_behav_eeg(:, ii) = temp_corr;    
        if (rem(ii,10) ==0)
            disp(['time bin ' num2str(ii) ' processed'])
        end 
end
%% setting variables for plots
variables = {'Liking', 'Valence', 'Arousal', 'Complexity', 'Familiarity', 'Style', 'Category'};
feature_colors = containers.Map(variables, {...
    [0.8500, 0.3250, 0.0980], ...  
    [0, 0.4470, 0.7410], ...       
    [0.4660, 0.6740, 0.1880], ...  
    [0.9290, 0.6940, 0.1250], ...  
    [0.4940, 0.1840, 0.5560], ...  
    [0.3010, 0.7450, 0.9330], ...  
    [0.6350, 0.0780, 0.1840] ...   
});
lineStyles = {'-', '-.'};
lineWidths = [2, 1, 3];

%% 
%%============= ADD STYLE AND CATEGORY PREDICTORS (BINARY) ===============

% STYLE predictor: 0 if same style, 1 if different style (column 3 in data_VAPS)
style_pred = [];
for ii = 1:999 
    st_1 =  data_VAPS.data(ii,3);
   for  jj = 1:999 
    st_2 =  data_VAPS.data(jj,3);
    if st_1 == st_2
        style_pred(ii,jj) = 0; 
    else 
        style_pred(ii,jj) = 1;
   end 
   end 
end 
style_predictor = squareform(style_pred)';
%% 
% CATEGORY predictor: 0 if same category, 1 if different category (column 4)
cat_pred = [];
for ii = 1:999 
    ct_1 =  data_VAPS.data(ii,4);
   for  jj = 1:999 
    ct_2 =  data_VAPS.data(jj,4);
    if ct_1 == ct_2
        cat_pred(ii,jj) = 0;
    else 
        cat_pred(ii,jj) = 1;
   end 
   end 
end 
cat_predictor = squareform(cat_pred)';
%% 
for kk = 1:2
    for ii =1: size(ds_corr,3)
    eeg_pred = squareform(squeeze(ds_corr(:,:,ii)))'; 
    if kk == 1
        corr_eeg(1,ii)= atanh(cosmo_corr(style_predictor,eeg_pred,'Spearman')); 
    else 
        corr_eeg(2,ii)= atanh(cosmo_corr(cat_predictor,eeg_pred,'Spearman'));
    end 
        if (rem(ii,50) ==0)
            disp(['time bin ' num2str(ii) ' processed'])
        end 
    end 
end 
%%
% Append style & category rows below the 5 behavioral dims
corr_behav_eeg_final = [corr_behav_eeg ; corr_eeg];
legend_handles = [];
%%===================== PLOT: SIMPLE CORRELATIONS =========================
h = plot(data_subset_bs.time, corr_behav_eeg_final(1,:),...
    "Color", feature_colors("Liking"), "LineStyle","-", ...
    "LineWidth", lineWidths(2), "DisplayName", "Liking");
legend_handles = [legend_handles, h];
hold on 
h = plot(data_subset_bs.time, corr_behav_eeg_final(2,:),...
    "Color", feature_colors("Valence"), "LineStyle","-", ...
    "LineWidth", lineWidths(2), "DisplayName", "Valence");
legend_handles = [legend_handles, h];
hold on 
h = plot(data_subset_bs.time, corr_behav_eeg_final(3,:),...
    "Color", feature_colors("Arousal"), "LineStyle","-", ...
    "LineWidth", lineWidths(2), "DisplayName", "Arousal");
legend_handles = [legend_handles, h];
hold on 
h = plot(data_subset_bs.time, corr_behav_eeg_final(4,:),...
    "Color", feature_colors("Complexity"), "LineStyle","-", ...
    "LineWidth", lineWidths(2), "DisplayName", "Complexity");
legend_handles = [legend_handles, h];
hold on 
h = plot(data_subset_bs.time, corr_behav_eeg_final(5,:),...
    "Color", feature_colors("Familiarity"), "LineStyle","-", ...
    "LineWidth", lineWidths(2), "DisplayName", "Familiarity");
legend_handles = [legend_handles, h];
hold on 
h = plot(data_subset_bs.time, corr_behav_eeg_final(6,:),...
    "Color", feature_colors("Style"), "LineStyle","-", ...
    "LineWidth", lineWidths(2), "DisplayName", "Style");
legend_handles = [legend_handles, h];
hold on 
h = plot(data_subset_bs.time, corr_behav_eeg_final(7,:),...
    "Color", feature_colors("Category"), "LineStyle","-", ...
    "LineWidth", lineWidths(2), "DisplayName", "Category");
legend_handles = [legend_handles, h];
hold on 
xline(0,'--'); hold on ; yline(0,'--');
xlabel('Time');
ylabel('Correlation');
legend(legend_handles, 'Location', 'best');
grid off;
ylim([min(corr_behav_eeg_final(:)), max(corr_behav_eeg_final(:))]);
xlim([min(data_subset_bs.time), max(data_subset_bs.time)]);
title(strcat(string(participant_num), 'th participant - Correlation'));
% Save plot
plotpath = fullfile(foldername, strcat('simple_corr_eeg_', string(participant_num), '_graph'));
savefig(plotpath);
%%
% Combine all 7 predictors: [5 behavioral dims, style, category]
feature_matrix = [pred_behav, style_predictor, cat_predictor];
num_vars = 7;
%%================ Partial Correlation =========================
partial_corr_behav_eeg = zeros(num_vars, numTimeBins);

for time_bin = 1:numTimeBins
    eeg_pred = squareform(squeeze(ds_corr(:,:, time_bin)))';
    
    for feature_idx = 1: num_vars
        current_feature = feature_matrix(:, feature_idx);
        all_feature_matrix = setdiff(1:num_vars, feature_idx);
        partial_corr_behav_eeg(feature_idx, time_bin) = partialcorr(current_feature, eeg_pred, feature_matrix(:, all_feature_matrix));
    end 
    if (rem(time_bin,10) ==0)
            disp(['time bin ' num2str(time_bin) ' processed'])
    end 
end
%%
figure;
hold on;
legend_handles = [];
for feature_idx = 1 : num_vars
    h = plot(data_subset_bs.time,partial_corr_behav_eeg(feature_idx, :), ...
        'Color', feature_colors(variables{feature_idx}), ...
        'DisplayName', variables{feature_idx});
    legend_handles = [legend_handles, h];
end

hold off;

xline(0, '--'); 
yline(0, '--');
ylim([min(partial_corr_behav_eeg(:)), max(partial_corr_behav_eeg(:))]);
xlim([min(data_subset_bs.time), max(data_subset_bs.time)]);
xlabel('Time');
ylabel('Partial Correlation');
title(strcat(string(participant_num), 'th participant - Partial Correlation'));
legend(legend_handles, 'Location', 'best');
grid off;
plotpath = fullfile(foldername, strcat('partial_corr_behav_eeg_', string(participant_num), '_graph'));
savefig(plotpath);

%%
% Save mat file for partial correlation
partial_file_name = fullfile(foldername, strcat('partial_corr_behav_eeg_', string(participant_num)));
save(partial_file_name, 'partial_corr_behav_eeg');

%%
% save mat file for simple correlation
simple_file_name = fullfile(foldername, strcat('simple_corr_eeg_', string(participant_num)));
save(simple_file_name, 'corr_behav_eeg_final');
%%
%%=============== COMPUTE STYLE‐SPECIFIC AND CATEGORY‐SPECIFIC CORRELATIONS ===============
% This section computes, for each style/category label, the Spearman‐based 
% “liking” representational dissimilarity (RDM) vs. EEG RDM correlation 
% (Fisher‐z transformed) across all time bins. 

styles = unique(data_VAPS.data(:, 3));
categories = unique(data_VAPS.data(:, 4));
num_styles = length(styles);
num_category = length(categories);
num_time_bins = size(ds_corr, 3);
%%
% Preallocate result matrices
corr_style_liking = zeros(num_styles, num_time_bins);
corr_category_liking = zeros(num_category, num_time_bins);
%%
%%-------------------------- STYLE‐SPECIFIC CORRELATIONS --------------------------

for s = 1: num_styles
    % Identify which stimuli belong to style #s 
    style_idx = data_VAPS.data(:, 3) == styles(s);
    % Extract the “liking” difference submatrix for this style
    diff_pred_styles = diff_pred(style_idx, style_idx);
    liking_values = squareform(diff_pred_styles)';
    
    % For each time bin, compute Pearson correlation between liking_values and EEG RDM
    for t = 1: num_time_bins
        eeg_pred = squareform(squeeze(ds_corr(style_idx, style_idx, t)))';
        corr_style_liking(s, t) = atanh(cosmo_corr(liking_values, eeg_pred, 'Pearson'));
    end
    disp(['Style ', num2str(styles(s)), ' processed']);
end

%%
%%------------------------ CATEGORY‐SPECIFIC CORRELATIONS ------------------------

for c = 1: num_category
    % Identify which stimuli belong to category #c
    cat_idx = data_VAPS.data(:, 4) == categories(c);
     % Extract “liking” difference submatrix for this category
    diff_pred_cat = diff_pred(cat_idx, cat_idx);
    liking_values = squareform(diff_pred_cat)';

    % For each time bin, compute Pearson correlation between liking_values and EEG RDM
    for t = 1: num_time_bins
        eeg_pred = squareform(squeeze(ds_corr(cat_idx, cat_idx, t)))';
        corr_category_liking(c, t) = atanh(cosmo_corr(liking_values, eeg_pred, 'Pearson'));
    end
    disp(['Category ', num2str(categories(c)), ' processed']);
end
%%
% initialize variable for plotting graph
style_name = ["Renaissance and Mannerism", "Baroque and Rococo", "Idealistic tendencies", ...
              "Realistic tendencies I", "Impressionistic tendencies", "Postimpressionistic tendencies", ...
              "Expressionistic tendencies", "Cubistic tendencies", "Realistic tendencies II", ...
              "Surrealistic tendencies", "Constructivist tendencies", "Abstract expressionistic tendencies", ...
              "Informal tendencies"];

category_name = ["Scenes", "Portrait", "Landscape", "Still life", "Toward Abstraction"];
colors_styles = lines(num_styles);
colors_categories = lines(num_category);
%%

%%================== PLOT: STYLE‐BY‐STYLE TIMECOURSES =======================
num_rows = ceil(sqrt(num_styles));
num_cols = ceil(num_styles / num_rows);

figure;
for i = 1:num_styles
    subplot(num_rows, num_cols, i);
    h = plot(data_subset_bs.time, corr_style_liking(i, :),...
        'Color', colors_styles(i, :), 'LineWidth', 1.5);
    hold on;
    mean_corr = mean(corr_style_liking(i, :)); 
    plot([min(data_subset_bs.time), max(data_subset_bs.time)],...
        [mean_corr, mean_corr], 'k--', 'LineWidth', 1.5,...
        'DisplayName', 'Mean Correlation');
    
    xline(0, '--'); 
    yline(0, '--');
    ylim([min(corr_style_liking(:)), max(corr_style_liking(:))]);
    xlim([min(data_subset_bs.time), max(data_subset_bs.time)]);
    
    title(style_name(i), 'FontSize', 10);
    xlabel('Time');
    ylabel('Correlation');
    grid off;
end

% save plot
sgtitle('Liking Correlation of Different Styles Over Time');
plotpath = fullfile(foldername, strcat('eeg_aesthetics_', string(participant), '\20Hz\corr_style_liking_', string(participant), '_graph.png'));
set(gcf, 'Position', [100, 100, 1200, 600]); 
exportgraphics(gcf, plotpath, 'Resolution', 300); 
%%
%%================ PLOT: CATEGORY‐BY‐CATEGORY TIMECOURSES ==================
num_rows = ceil(sqrt(num_category)); 
num_cols = ceil(num_category / num_rows);

figure;
for i = 1:num_category
    subplot(num_rows, num_cols, i);
    plot(data_subset_bs.time, corr_category_liking(i, :),...
        'Color', colors_categories(i, :), 'LineWidth', 1.5);
    hold on;
    mean_corr = mean(corr_category_liking(i, :)); 
    plot([min(data_subset_bs.time), max(data_subset_bs.time)],...
        [mean_corr, mean_corr], 'k--', 'LineWidth', 1.5,...
        'DisplayName', 'Mean Correlation');
    
    xline(0, '--'); 
    yline(0, '--');
    ylim([min(corr_category_liking(:)), max(corr_category_liking(:))]);
    xlim([min(data_subset_bs.time), max(data_subset_bs.time)]);
    
    title(category_name(i), 'FontSize', 10); 
    xlabel('Time');
    ylabel('Correlation');    
    grid off;
end

sgtitle('Liking Correlation of Different Categories Over Time');
plotpath = fullfile(foldername, strcat('eeg_aesthetics_', string(participant), '\20Hz\corr_category_liking_', string(participant), '_graph.png'));
set(gcf, 'Position', [100, 100, 1200, 600]); 
exportgraphics(gcf, plotpath, 'Resolution', 300); 
%%
%%----------------- SAVE CORRELATION MATRICES TO .MAT FILES ------------------
corr_style_liking_filename = fullfile(foldername, strcat('eeg_aesthetics_', string(participant), '\20Hz\corr_style_liking_', string(participant)));
save(corr_style_liking_filename, 'corr_style_liking');

corr_category_liking_filename = fullfile(foldername, strcat('eeg_aesthetics_', string(participant), '\20Hz\corr_category_liking_', string(participant)));
save(corr_category_liking_filename, 'corr_category_liking');