%%

clc
close all
clear all
 addpath(genpath('/Users/ayaygoya/Documents/VaishInternship/CoSMoMVPA/'));

logfile = load('/Users/ayaygoya/Documents/VaishInternship/Data/log/RSVP_eeg_s1.mat');
% load vaps data

data_VAPS = load('/Users/ayaygoya/Documents/VaishInternship/database_VAPS.mat');
 %% 
for ii = 1: length(logfile.dat.new_stim)
    
    if logfile.dat.trigs(ii) ==1
        image_presented = logfile.dat.new_stim{ii};
        
        %inputString = 'image_12.jpg';
        
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
        %disp(['Extracted number: ', num2str(extractedNumber)]);
        
        trials_data(ii,:) = data_VAPS.data(extracted_number,:);
    end
end

%%
ft_defaults()

fileName = ('/Users/ayaygoya/Documents/VaishInternship/Data/eeg/eeg_aesthetics_0001.eeg');
%Define Events
cfg=[]; 
cfg.dataset=fileName;
%cfg.trialdef.eventtype='Stimulus';
cfg.trialdef.eventtype= 'Stimulus';

%cfg.trialdef.eventvalue = ?

cfg.trialdef.prestim=0.1;
cfg.trialdef.poststim=1;
cfg=ft_definetrial(cfg);

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
%%
cfg=[];
cfg.showlabel='yes';
cfg.method='summary'; 
cfg.keepchannel='no';
data_clean=ft_rejectvisual(cfg,data);
%%
% downsampling data to 100Hz
cfg = []
cfg.resamplefs = 100; 
data_selected = ft_resampledata(cfg,data_clean)
%%
% Run ICA
cfg = [];
cfg.method = 'runica';
comp = ft_componentanalysis(cfg, data_selected);
%%
% Inspect ICA components
cfg = [];
cfg.component = 1:20; 
cfg.layout = 'easycapM1.mat'; 
ft_topoplotIC(cfg, comp);

cfg = [];
cfg.viewmode = 'component';
ft_databrowser(cfg, comp);
%%
% Remove blink components
cfg = [];
cfg.component = [1 2]; % components to be removed
data_clean = ft_rejectcomponent(cfg, comp, data_selected);
%%
cfg= [];
cfg.trials = find(data_clean.trialinfo(:,1) ==1);
data_subset = ft_selectdata(cfg,data_clean);

cfg=[];
cfg.keeptrials = 'yes';
data_subset = ft_timelockanalysis(cfg,data_subset);
%%
% data_classify = data_subset.trial;
% 
% 
% cfg=[]
% cfg.classifier = 'multiclass_lda';
% cfg.metric = 'accuracy';
% cfg.cv = 'kfold';
% cfg.k = 5;
% cfg.repeat = 5;
% cfg.preprocess = {'undersample'};
% cfg.undersample_test_set = 1
% [auc,results] = mv_classify_across_time(cfg, data_classify , data_subset.trialinfo(:,5))
% 
% 
% figure
% plot(data_subset.time, auc)
% hold on
% xline(0,'--');
% hold on
% yline(0.20,'--')
% xlabel('time')
% ylabel('Decoding Accuracy')
% title('Decoding Art category: 5 class')
% 
% 
% %%
% 
% cfg=[]
% cfg.classifier = 'multiclass_lda';
% cfg.metric = 'accuracy';
% cfg.cv = 'kfold';
% cfg.k = 5;
% cfg.repeat = 1;
% cfg.preprocess = {'undersample'};
% cfg.undersample_test_set = 1
% [auc,results] = mv_classify_across_time(cfg, data_classify , data_subset.trialinfo(:,4))
% 
% 
% figure
% plot(data_subset.time, auc)
% hold on
% xline(0,'--');
% hold on
% yline(0.076,'--')
% xlabel('time')
% ylabel('Decoding Accuracy')
% title('Decoding Artistic style: 13 classes')
% %% 
% cfg=[]
% cfg.trials = find(data_subset.trialinfo(:,8)==1 | data_subset.trialinfo(:,8)==2)
% data_ort = ft_selectdata(cfg,data_subset)
% 
% cfg=[]
% cfg.classifier = 'lda';
% cfg.metric = 'accuracy';
% cfg.cv = 'kfold';
% cfg.k = 5;
% cfg.repeat = 1;
% cfg.preprocess = {'undersample'};
% cfg.undersample_test_set = 1
% [auc,results] = mv_classify_across_time(cfg, data_ort.trial , data_ort.trialinfo(:,8))
% 
% %%
% figure
% plot(data_subset.time, auc)
% hold on
% xline(0,'--');
% hold on
% yline(0.50,'--')
% xlabel('time')
% ylabel('Decoding Accuracy')
% title('Decoding Art orientation: 2 classes')


%% 
% baseline correction 
cfg=[]
cfg.baseline = [-0.05 0]
data_subset_bs = ft_timelockbaseline(cfg,data_subset)

cfg = []
cfg.resamplefs = 100; 
data_subset_bs = ft_resampledata(cfg,data_subset_bs)

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
    
    %ds_avg_samples.samples = double(ds_avg_samples);
    
    ds_corr(:,:,n_bin) = 1 - cosmo_corr(ds_avg_samples.samples','Spearman');
    if (rem(n_bin,50) ==0)
        disp(['time bin ' num2str(n_bin) ' processed'])
    end 
end 


%% 

cond_colms = [8, 10, 12, 14,16];
for n_cond = 1:5
for ii = 1:999 
    
    %  for liking 8 for valence 10 %  arousal 12, complexity, 14, Familarity
    % 16
    
    
    pred_val_1 = data_VAPS.data(ii,cond_colms(n_cond));
    
    for jj = 1:999
        
        pred_val_2 = data_VAPS.data(jj,cond_colms(n_cond));
        
        diff_pred(ii,jj) = abs(pred_val_1 - pred_val_2);
        
    end 
    
    
end

pred_behav(:,n_cond) = squareform(diff_pred)';

end 

%% 
for kk = 1:5
    for ii =1: size(ds_corr,3)
    
    eeg_pred = squareform(squeeze(ds_corr(:,:,ii)))'; 
    
    corr_behav_eeg(kk,ii)= atanh(corr(pred_behav(:,kk),eeg_pred,'type','Spearman')); 
    
        if (rem(ii,50) ==0)
            disp(['time bin ' num2str(ii) ' processed'])
        end 

    end 

end 
% for ii =1: size(ds_corr,3)
%     
%     eeg_pred = squareform(squeeze(ds_corr(:,:,ii)))'; 
%     
%     corr_behav_eeg_2(1,ii)= atanh(corr(eeg_pred,pred_behav,'type','Spearman')); 
%     
%     if (rem(ii,50) ==0)
%         disp(['time bin ' num2str(ii) ' processed'])
%     end 
%     
% end 

% 
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
title('Second Participant')

%% 

cfg = []
cfg.layout = 'easycapM1.lay'
ft_multiplotER(cfg,data_subset)

%% 

style_pred = []
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
cat_pred = []
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
    if kk == 1;
        corr_eeg(1,ii)= atanh(corr(style_predictor,eeg_pred,'type','Spearman')); 
    else 
        corr_eeg(2,ii)= atanh(corr(cat_predictor,eeg_pred,'type','Spearman')); 

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
title('Second Participant')
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
title('Second Participant')
