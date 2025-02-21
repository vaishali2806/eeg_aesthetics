%% Define some stuff!
tic
clc %Remove old stuff
clear all
close all
load('C:\Users\JLU-SU\Desktop\Stimilu task\data-ort.mat');
try

    start_time=GetSecs; %Start timer

    addpath('C:\Users\JLU-SU\Desktop\Stimilu task\art_all_renamed\'); %Add the path to stimuli

    cfg.eeg_mode=0; %toggle EEG mode (set to 0 to run it without the EEG)

    % Added for the trigger
    if cfg.eeg_mode==1
        SerialPortObj=serial('COM4', 'TimeOut',1);%COM4 is the virtual serial output in the device manager
        SerialPortObj.BytesAvailableFcnMode='byte';
        SerialPortObj.BytesAvailableFcnCount=1;
        fopen(SerialPortObj);
        fwrite(SerialPortObj, 0, 'sync');
    end

    cfg.windowcolor=[125,125,125]; %Colors and Fonts
    cfg.textcolor=[255,255,255];
    cfg.fixcolor=[225,25,225];
    cfg.fontname='Arial';
    cfg.fontsize=20;

    cfg.nStim=999;
    cfg.nRep=12;
    cfg.trialamount=cfg.nRep*cfg.nStim;  % n Stimuli, n Repetitions
    cfg.blocksize = 200;

    cfg.picsize_hor=300; %Size of the images (this needs to be square for the ratings, pix size is fixed below)
    cfg.picsize_vert=300;

    cfg.deg_object=10; %desired degree (visual angle) of the object

    [w,h]=Screen('DisplaySize',0); %Get display size in mm
    cfg.screen_hor=w/10; %Convert to cm (horizontal)
    cfg.screen_vert=h/10; %Convert to cm (vertial)

    cfg.screen_dist=60; %Distance to screen in cm

    w=Screen('Resolution',0); %Get screen properties
    cfg.screensize=[w.width,w.height]; %./2; %Store horizontal/vertical pixels (i.e., resolution)

    cfg.time_stim=.05; %Duration of the Stimulus
    cfg.time_wait=1; %quick wait between after the stimulus
    % cfg.time_iti=0.15; %Duration of the inter-trial interval (will be jittered by +/-200ms)

    KbName('UnifyKeyNames'); %Fix Mac/Windows differences (hopefully)
    cfg.escapeKey=KbName('q'); %Get key Codes (indices in the KeyCode Matrix) for a few keys
    cfg.spaceKey=KbName('space');

    dat.subjcode=input('Enter subject number: '); %Ask for subject number
    dat.sub_gender = input('Enter subject gender (1: Male 2: Female 3:other ) : ');
    dat.sub_age = input('Enter subject age : ');
    rng(dat.subjcode); %Set random numbers to be truly random (this can throw an error -> remove it then, it does the same more clumsily below)

    %% Screen Setup
    for i=1:dat.subjcode  %This is a workaround if the random number generator doesn't work (e.g., for old matlab versions) to make the random sequences different for all people
        rand(1,1);randperm(10);Shuffle(1:10);randi(10);
    end

    Screen('Preference', 'SkipSyncTests', 0);
    [mainwindow,screen_rect]=Screen('OpenWindow',0,cfg.windowcolor);

    newTextSize=Screen('TextSize', mainwindow, [cfg.fontsize]);  %Specify Font Size
    Screen('TextFont', mainwindow, [cfg.fontname]);  %Specify Font

    cfg.monitorFrameRate=FrameRate(mainwindow); %Get the frame rate of the screen
    cfg.monitorFlipInterval=Screen('GetFlipInterval',mainwindow); %Get the interval between flips (is a function of the framerate, e.g., 16ms for 60hz)

    real_stim=1000*cfg.time_stim; %For the stimulus
    stim_frames=real_stim/(1000/cfg.monitorFrameRate);
    time_stim=cfg.monitorFlipInterval*(stim_frames-0.5);

    % real_iti=1000*cfg.time_iti; %For the stimulus
    % iti_frames=real_iti/(1000/cfg.monitorFrameRate);
    % time_iti=cfg.monitorFlipInterval*(iti_frames-0.5);

    clear hidecursor;
    %HideCursor;

    screensize= cfg.screensize;
    center = cfg.screensize./2;  %Determine center coordinates
    center_hor=center(1);
    center_vert=center(2);

    % This is for the stimulus size:
    picsizecm=tan((cfg.deg_object/2)*2*pi/360)*cfg.screen_dist;  %Centimeters size of the picture
    picsize=[picsizecm*cfg.screensize(1)/cfg.screen_hor,picsizecm*cfg.screensize(1)/cfg.screen_hor*(cfg.picsize_vert/cfg.picsize_hor)];  %Pixels of the picture

    stimpositionmid=center; %-[0,screensize(1)*0.05];
    fixposition=stimpositionmid;

    %all-purpose quadratic rectangle
    picrectangle1=[stimpositionmid(1)-picsize(1),stimpositionmid(2)-picsize(2)];
    picrectangle2=[stimpositionmid(1)+picsize(1),stimpositionmid(2)+picsize(2)];
    picrectangle=[picrectangle1,picrectangle2];  %Rectangle for the Picture

    %rectangle for stimulus
    aspect=1;
    picrectangle1=[stimpositionmid(1)-aspect*picsize(1),stimpositionmid(2)-picsize(2)];
    picrectangle2=[stimpositionmid(1)+aspect*picsize(1),stimpositionmid(2)+picsize(2)];
    stimrectangle=[picrectangle1,picrectangle2];  %Rectangle for the Picture

    %rectangle for rating wheel
    sf=1.35;
    picrectangle1=[stimpositionmid(1)-picsize(1)*sf,stimpositionmid(2)-picsize(2)*sf];
    picrectangle2=[stimpositionmid(1)+picsize(1)*sf,stimpositionmid(2)+picsize(2)*sf];
    wheelrectangle=[picrectangle1,picrectangle2];  %Rectangle for the heel

    stimpositionmid=center;

    %rectangles for rating bars
    bar1=[stimpositionmid(1)-screensize(1)/3,stimpositionmid(2)+screensize(2)/4];
    bar2=[stimpositionmid(1)+screensize(1)/3,stimpositionmid(2)+screensize(2)/4+screensize(2)/20];
    bar=[bar1,bar2];
    for i=1:7
        subbar1=[stimpositionmid(1)-screensize(1)/3+(bar2(1)-bar1(1))/7*(i-1),stimpositionmid(2)+screensize(2)/4];
        subbar2=[stimpositionmid(1)-screensize(1)/3+(bar2(1)-bar1(1))/7*i,stimpositionmid(2)+screensize(2)/4+screensize(2)/20];
        subbar{i}=[subbar1,subbar2];
    end
    for i=1:7
        subbar1=[stimpositionmid(1)-screensize(1)/3+(bar2(1)-bar1(1))/7*(i-1),stimpositionmid(2)+screensize(2)/4+1.5*screensize(2)/20];
        subbar2=[stimpositionmid(1)-screensize(1)/3+(bar2(1)-bar1(1))/7*i,stimpositionmid(2)+screensize(2)/4+2.5*screensize(2)/20];
        subbar_low{i}=[subbar1,subbar2];
    end
    for i=1:7
        subbar1=[stimpositionmid(1)-screensize(1)/3+(bar2(1)-bar1(1))/7*(i-1),stimpositionmid(2)+screensize(2)/4+3*screensize(2)/20];
        subbar2=[stimpositionmid(1)-screensize(1)/3+(bar2(1)-bar1(1))/7*i,stimpositionmid(2)+screensize(2)/4+4*screensize(2)/20];
        subbar_lowlow{i}=[subbar1,subbar2];
    end

    %% Stimulus Randomisation
    num_blocks = round(cfg.nStim * cfg.nRep / cfg.blocksize);

    t1=[];
    for b=1:cfg.nRep
        t1=[t1;randperm(cfg.nStim)'];
    end

    for trial=1:cfg.trialamount
        %define stimulus name
        stim_nr=t1(trial);
        dat.ort{trial} = data_ort(stim_nr,2);
        dat.stim{trial}=['image_',num2str(stim_nr),'.jpg'];
        dat.aspect_ratio{trial} = data_ort(stim_nr,5);
    end

    dat = add_distractors_faster(dat, cfg.blocksize);

    dat = generate_mask_indices(dat,200);

    dat = generate_trigs(dat);

    % Pre-load all images
    timestart = GetSecs;
    patterns = {'dist_', 'image_', 'Masks_'};
    imageFiles = dir('C:\Users\JLU-SU\Desktop\Stimilu task\art_all_renamed\*.jpg');
    imageData = struct();
    fprintf('Loading images...\n');
    for i = 1:numel(imageFiles)
        progress = i / length(imageFiles);
    
        % Print the progress to the command window
        fprintf('Loading: %.2f%%\n', progress * 100);
        fileName = imageFiles(i).name;
        
        % Extract image number using regex
        imgNum = regexp(fileName, '\d+', 'match', 'once');
        
        if isempty(imgNum)
            continue; % Skip files without a number
        end
    
        imgNum = str2double(imgNum);
        
        % Determine the type based on the filename prefix
        for p = 1:numel(patterns)
            if contains(fileName, patterns{p})
                imageType = patterns{p}(1:end-1); % Remove underscore
                break;
            end
        end
    
        % Store the image in a structured format
        if ~isfield(imageData, imageType)
            imageData.(imageType) = containers.Map('KeyType', 'double', 'ValueType', 'any');
        end
        imageData.(imageType)(imgNum) = imread(fileName);
    end
    
    fprintf('All images loaded successfully!\n');
    timeend = GetSecs;
    image_load_time = timeend-timestart;
    fprintf('Elapsed time: %.2f%\n', timeend-timestart);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Start Presentation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    DrawFormattedText(mainwindow,['instructions \n \n',...
        'please keep central fixation \n \n ', ...
        'then use the mouse to answer the questions. \n \n press space to proceed...'],'center','center',cfg.textcolor);
    Screen('Flip',mainwindow);
    KbWait(-1);

    trial_timing = [];
    for trial=1: length(dat.new_stim)%cfg.trialamount

        escapeKey = KbName('escape');
        [ keyIsDown, seconds, keyCode ] = KbCheck(-1);
        if keyCode(escapeKey) == 1
            break;

            sca;
        end

        %stimulus on
        stimFilename = dat.new_stim{trial}; % e.g., 'image_1.jpg'
        stimNum = str2double(regexp(stimFilename, '\d+', 'match', 'once'));

        if contains(stimFilename, 'dist_')
            stimType = 'dist';
        elseif contains(stimFilename, 'image_')
            stimType = 'image';
        else
            error('Unknown stimulus type: %s', stimFilename);
        end
        % Check if the image type exists in `imageData`
        if ~isfield(imageData, stimType)
            error('Missing image category: %s', stimType);
        end
        % Check if the extracted number exists in the `containers.Map`
        if ~imageData.(stimType).isKey(stimNum)
            error('Image %s_%d not found in preloaded imageData.', stimType, stimNum);
        end
        img = imageData.(stimType)(stimNum);       % Preloaded image
        
        % Get the mask number from dat.mask_indices
        maskNum = dat.mask_indices(trial);
        % Check if mask exists before accessing it
        if ~isfield(imageData, 'Masks') || ~imageData.Masks.isKey(maskNum)
            error('Mask image Masks_%d.jpg not found in preloaded imageData.', maskNum);
        end
        
        % Retrieve mask directly from preloaded imageData
        img_rand_orig = imageData.Masks(maskNum);  % Preloaded 'Masks' image
        
        % Apply desaturation function
        img_rand = desaturate_image1(img_rand_orig);
        
        % Convert images to Psychtoolbox textures
        stim_tex_rand = Screen('MakeTexture', mainwindow, img_rand);
        stim_tex = Screen('MakeTexture', mainwindow, img);
        
        % Draw to screen
        Screen('DrawTexture', mainwindow, stim_tex_rand, [], stimrectangle);

        if dat.new_ort{trial} == 1
            rect_size = [(stimrectangle(3) - stimrectangle(1)) (stimrectangle(4) - stimrectangle(2)) ];
            img_size = rect_size(2) * dat.new_asp{trial};
            diff = rect_size(1) - img_size;
            Screen('DrawTexture',mainwindow,stim_tex,[],[(stimrectangle(1) + diff/2)  stimrectangle(2) (stimrectangle(1)+ diff/2 + img_size)  stimrectangle(4)]);
        else
            rect_size = [(stimrectangle(3) - stimrectangle(1)) (stimrectangle(4) - stimrectangle(2)) ];
            img_size = rect_size(1) / dat.new_asp{trial};
            diff = rect_size(2) - img_size;
            Screen('DrawTexture',mainwindow,stim_tex,[],[stimrectangle(1)   (stimrectangle(2)+ diff/2) stimrectangle(3)  (stimrectangle(2) + img_size + diff/2)]);
        end
        Screen('DrawDots',mainwindow,fixposition,5,cfg.fixcolor);
        if trial == 1
            ON =Screen('Flip',mainwindow);
        end

          % Give a trigger once the stimuli is ON
              if cfg.eeg_mode == 1
                 fwrite(SerialPortObj,dat.trigs(trial),'sync');
                 pause(0.01);
                 fwrite(SerialPortObj,0,'sync');
              end

        flip_1_timing = ON;

        if trial ~= 1
            flip_duration(trial,2) = flip_1_timing - flip_2_timing;
        end

        Screen('DrawDots',mainwindow,fixposition,5,cfg.fixcolor);
        ON=Screen('Flip',mainwindow,ON + time_stim);

        flip_2_timing =ON;
        flip_duration(trial,1) = flip_2_timing - flip_1_timing;
        
        %

        if ~(rem(trial, cfg.blocksize))
            ON=Screen('Flip',mainwindow,ON);
            DrawFormattedText(mainwindow,['Instructions \n \n',...
                'please keep central fixation \n \n ', ...
                'press space to proceed...'],'center','center',cfg.textcolor);
            Screen('Flip',mainwindow);
            KbWait(-1);
        end
        HideCursor;
        %close texture
        Screen('Close',stim_tex_rand);
        Screen('Close',stim_tex);
    end
    sca
    if cfg.eeg_mode==1
        fclose(SerialPortObj);
        delete(SerialPortObj);
        clear SerialPortObj;
    end
    save(['RSVP_eeg_s',num2str(dat.subjcode)],'cfg','dat');
catch
    sca
    lasterr
end
toc
