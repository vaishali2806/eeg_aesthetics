%%

clc
close all
clear all
 addpath(genpath('D:\Internship\CoSMoMVPA-master'));

logfile = load('D:\Internship\Data-VG\log\RSVP_eeg_s3.mat');
% load vaps data
%%
data_VAPS = load('D:\Internship\database_VAPS.mat');
 %% 
 save_idx = 1;
for ii = 1: length(logfile.dat.new_stim)
    
    if logfile.dat.trigs(ii) ==1
        image_presented = logfile.dat.new_stim{ii};
        
        % Use a regular expression to extract the number
        pattern = '\d+';
        number = regexp(image_presented, pattern, 'match');
        
        % Convert the extracted number from cell to double
        if ~isempty(number)
            extracted_number = str2double(number{1});
        else
            extracted_number = [];
        end
        
        % Display the result        
        trials_data(save_idx,:) = data_VAPS.data(extracted_number,:);
        save_idx = save_idx + 1;
    end
end

%%
ft_defaults()

fileName = ('D:\Internship\Data-VG\eeg\eeg_aesthetics_0003.eeg');
%Define Events
cfg=[]; 
cfg.dataset=fileName;
cfg.trialdef.eventtype= 'Stimulus';

cfg.trialdef.prestim=0.1;
cfg.trialdef.poststim=1;
cfg=ft_definetrial(cfg);
%%
% Load Data
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
data.trialinfo  = [data.trialinfo trials_data];
%clearvars extracted_number ii image_presented number pattern trials_data logfile
%%
% downsampling data to 100Hz
cfg = [];
cfg.resamplefs = 100; 
data_selected = ft_resampledata(cfg,data);
%%
zero_var_channels = [];
for i = 1:length(data_selected.label)
    channel_variances = zeros(1, length(data_selected.trial));
    for j = 1:length(data_selected.trial)
        channel_data = data_selected.trial{j}(i, :);
        channel_variances(j) = var(channel_data);
    end
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
cfg=[];
cfg.showlabel='yes';
cfg.method='summary'; 
cfg.keepchannel='no';
data_clean=ft_rejectvisual(cfg,clean_data);
%clearvars data_selected
%%
% Run FastICA
cfg = [];
cfg.method = 'runica';
cfg.numcomponent = 30;
comp = ft_componentanalysis(cfg, data_clean);
%%
% Inspect ICA components
cfg = [];
cfg.component = 1:30; 
cfg.layout = 'easycap-M1.txt'; 
cfg.marker = 'off';
ft_topoplotIC(cfg, comp);
%%
cfg = [];
cfg.viewmode = 'component';
cfg.channel = [6 26 28 30];
cfg.layout    = 'easycap-M1.txt';
ft_databrowser(cfg, comp);
% %%
% Remove blink components
 cfg = [];
 cfg.component = [6];
 data_clean = ft_rejectcomponent(cfg, comp, data_clean);
%%
cfg= [];
cfg.trials = find(data_clean.trialinfo(:,1) ==1);
data_subset = ft_selectdata(cfg,data_clean);

cfg=[];
cfg.keeptrials = 'yes';
data_subset = ft_timelockanalysis(cfg,data_subset);
%clearvars data_clean
%% 
% baseline correction 
cfg=[];
cfg.baseline = [-0.05 0];
data_subset_bs = ft_timelockbaseline(cfg,data_subset);
%%
ds=cosmo_meeg_dataset(data_subset_bs);
ds.sa.targets = data_subset.trialinfo(:,2);
ds.sa.chunks = ones(length(data_subset.trialinfo),1);
%% 
for n_bin = 1:size(data_subset_bs.trial,3)
    ds_t = cosmo_slice(ds,ismember(ds.fa.time,n_bin),2);
    
    ds_pca = cosmo_map_pca(ds_t, 'pca_explained_ratio', 0.99, 'max_feature_count', 100000);
    ds_pca.sa.targets = data_subset.trialinfo(:,2);
    ds_pca.sa.chunks = ones(length(data_subset.trialinfo),1);
    ds_avg_samples = cosmo_average_samples(ds_pca);
    
    ds_corr(:,:,n_bin) = 1 - cosmo_corr(ds_avg_samples.samples','Spearman');
    if (rem(n_bin,50) ==0)
        disp(['time bin ' num2str(n_bin) ' processed'])
    end 
end

%% 
cond_colms = [8, 10, 12, 14,16];
for n_cond = 1:5
for ii = 1:999 
    %  for liking 8 for valence 10 %  arousal 12, complexity, 14, Familarity 16
    pred_val_1 = data_VAPS.data(ii,cond_colms(n_cond));
    for jj = 1:999
        pred_val_2 = data_VAPS.data(jj,cond_colms(n_cond));
        diff_pred(ii,jj) = abs(pred_val_1 - pred_val_2);
    end 
end
pred_behav(:,n_cond) = squareform(diff_pred)';
end 
%% 
numRowsPredBehav = size(pred_behav, 1);
numColsPredBehav = size(pred_behav, 2);
numTimeBins = size(ds_corr, 3);

corr_behav_eeg = zeros(numColsPredBehav, numTimeBins);

for ii =1: numTimeBins
    eeg_pred = squareform(squeeze(ds_corr(:,:,ii)))';
    temp_corr = zeros(numColsPredBehav, 1);
    for kk = 1:numColsPredBehav
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

cfg = [];
cfg.layout = 'easycapM1.lay';
ft_multiplotER(cfg,data_subset)

%% 
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
corr_behav_eeg_final = [corr_behav_eeg ; corr_eeg];
legend_handles = [];
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
ylabel('Simple Correlation');
legend(legend_handles, 'Location', 'best');
grid off;
ylim([min(corr_behav_eeg_final(:)), max(corr_behav_eeg_final(:))]);
xlim([min(data_subset_bs.time), max(data_subset_bs.time)]);
title('Thirteenth Participant - Simple Correlation')
%%
feature_matrix = [pred_behav, style_predictor, cat_predictor];
num_vars = 7;
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
title('Partial Correlation');
legend(legend_handles, 'Location', 'best');
grid off;

%%
% Computing correlations based on style and category 

% Extracting unique styles and categories
styles = unique(data_VAPS.data(:, 3));
categories = unique(data_VAPS.data(:, 4));

num_styles = length(styles);
num_category = length(categories);
num_time_bins = size(ds_corr, 3);
%%
corr_style_liking = zeros(num_styles, num_time_bins);
corr_category_liking = zeros(num_category, num_time_bins);
%%
% computing correlation for each style

for s = 1: num_styles
    style_idx = data_VAPS.data(:, 3) == styles(s);
    liking_values = data_VAPS.data(style_idx, 8);

    style_liking_pred = squareform(pdist(liking_values, 'cityblock'));

    for t = 1: num_time_bins
        eeg_pred = ds_corr(style_idx, style_idx, t);
        corr_style_liking(s, t) = atanh(cosmo_corr(style_liking_pred(:), eeg_pred(:), 'Spearman'));
    end
    disp(['Style ', num2str(styles(s)), ' processed']);
end

%%
% Computing Correlationg for each category 

for c = 1: num_category
    cat_idx = data_VAPS.data(:, 4) == categories(c);
    liking_values = data_VAPS.data(cat_idx, 8);

    cat_liking_pred = squareform(pdist(liking_values, 'cityblock'));

    for t = 1: num_time_bins
        eeg_pred = ds_corr(cat_idx, cat_idx, t);
        corr_category_liking(c, t) = atanh(cosmo_corr(cat_liking_pred(:), eeg_pred(:), 'Spearman'));
    end
    disp(['Category ', num2str(categories(c)), ' processed']);
end
%%
style_name = ["Renaissance and Mannerism", "Baroque and Rococo", "Idealistic tendencies", ...
              "Realistic tendencies I", "Impressionistic tendencies", "Postimpressionistic tendencies", ...
              "Expressionistic tendencies", "Cubistic tendencies", "Realistic tendencies II", ...
              "Surrealistic tendencies", "Constructivist tendencies", "Abstract expressionistic tendencies", ...
              "Informal tendencies"];

category_name = ["Scenes", "Portrait", "Landscape", "Still life", "Toward Abstraction"];
colors_styles = lines(num_styles);
colors_categories = lines(num_category);
%%

% figure; hold on;
% legend_styles = [];
% for i = 1:num_styles
%    h = plot(data_subset_bs.time, corr_style_liking(i, :),...
%        'Color', colors_styles(i, :),...
%        'LineWidth', 1.5,...
%        'DisplayName', style_name(i));
%    legend_styles = [legend_styles, h];
% end
% plot(data_subset_bs.time, mean(corr_style_liking, 1), 'k', 'LineWidth', 2, 'DisplayName', 'Mean Correlation');
% 
% xline(0, '--'); 
% yline(0, '--');
% ylim([min(corr_style_liking(:)), max(corr_style_liking(:))]);
% xlim([min(data_subset_bs.time), max(data_subset_bs.time)]);
% xlabel('Time');
% ylabel('Correlation');
% title('Correlation of Different Styles');
% legend(legend_styles, 'Location', 'best');
% grid off;
% hold off;
% %%
% figure;
% hold on;
% legend_cat = [];
% for i = 1:num_category
%    h = plot(data_subset_bs.time, corr_category_liking(i, :),...
%        'Color', colors_categories(i, :),...
%        'LineWidth', 1.5,...
%        'DisplayName', category_name(i));
%    legend_cat = [legend_cat, h];
% end
% 
% hold off;
% 
% xline(0, '--'); 
% yline(0, '--');
% ylim([min(corr_category_liking(:)), max(corr_category_liking(:))]);
% xlim([min(data_subset_bs.time), max(data_subset_bs.time)]);
% xlabel('Time');
% ylabel('Correlation');
% title('Correlation of Different Categories');
% legend(legend_cat, 'Location', 'best');
% grid off;
%%

% Style Graph
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

sgtitle('Liking Correlation of Different Styles Over Time');
%%

%Category Graph
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
