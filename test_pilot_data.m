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
        trials_data(ii,:) = data_VAPS.data(extracted_number,:);
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
%% 
plot(data_subset_bs.time, corr_behav_eeg(1,:))
hold on 
plot(data_subset_bs.time, corr_behav_eeg(2,:))
hold on 
plot(data_subset_bs.time, corr_behav_eeg(3,:))
hold on 
plot(data_subset_bs.time, corr_behav_eeg(4,:))
hold on 
plot(data_subset_bs.time, corr_behav_eeg(5,:))
hold on 
xline(0,'--'); hold on ; yline(0,'--')
xlabel('Time')
ylabel('Correlation')
legend('Liking','Valence','Arousal', 'Complexity','Familarity')
xlim([-0.05 1.0])
title('Third Participant')

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
plot(data_subset_bs.time, corr_eeg(1,:))
hold on 
plot(data_subset_bs.time, corr_eeg(2,:))
hold on 
xline(0,'--'); hold on ; yline(0,'--')
xlabel('Time')
ylabel('Correlation')
legend('Style', 'Category')
xlim([-0.05 1.0])
title('Third Participant')
%%
plot(data_subset_bs.time, corr_behav_eeg(1,:))
hold on 
plot(data_subset_bs.time, corr_behav_eeg(2,:))
hold on 
plot(data_subset_bs.time, corr_behav_eeg(3,:))
hold on 
plot(data_subset_bs.time, corr_behav_eeg(4,:))
hold on 
plot(data_subset_bs.time, corr_behav_eeg(5,:))
hold on 
plot(data_subset_bs.time, corr_eeg(1,:))
hold on 
plot(data_subset_bs.time, corr_eeg(2,:))
hold on 
xline(0,'--'); hold on ; yline(0,'--')
xlabel('Time')
ylabel('Correlation')
legend('Liking','Valence','Arousal', 'Complexity','Familarity', 'Style','Category')
xlim([-0.05 1.0])
title('Third Participant')
%%
variables = {'Liking', 'Valence', 'Arousal', 'Complexity', 'Familiarity', 'Style', 'Category'};
feature_matrix = [pred_behav, style_predictor, cat_predictor];
num_vars = 7;
partial_corr_behav_eeg = zeros(num_vars, num_vars, numTimeBins);

for feature_idx_1 = 1:num_vars
    for feature_idx_2 = 1:num_vars
        if feature_idx_2 ~= feature_idx_1
            for time_bin = 1: numTimeBins
                eeg_pred = squareform(squeeze(ds_corr(:,:, time_bin)))';
                current_feature = feature_matrix(:, feature_idx_1);
                control_vars = setdiff(1:num_vars, [feature_idx_1, feature_idx_2]);
                partial_corr_behav_eeg(feature_idx_1, feature_idx_2, time_bin) = ...
                    partialcorr(current_feature, eeg_pred, feature_matrix(:, control_vars));
            end
            if rem(time_bin, 10) == 0
                disp(['Time bin ' num2str(time_bin) ' processed for variable ' num2str(feature_idx_1) ' and variable ' num2str(feature_idx_2)]);
            end
        end
    end
end
%%
corr_behav_eeg_final = [corr_behav_eeg ; corr_eeg];
for main_var = 1:num_vars
    figure;
    plot(data_subset_bs.time, corr_behav_eeg_final(main_var, :), 'DisplayName', variables{main_var});
    hold on;
    for removed_var = 1:num_vars
        if main_var ~= removed_var
            plot(data_subset_bs.time, squeeze(partial_corr_behav_eeg(main_var, removed_var, :)), 'DisplayName', ['Removed ' variables{removed_var}]);
        end
    end
    xline(0, '--'); 
    yline(0, '--');
    xlabel('Time');
    ylabel('Partial Correlation');
    title(['Partial Correlation for ' variables{main_var}]);
    legend;
end

