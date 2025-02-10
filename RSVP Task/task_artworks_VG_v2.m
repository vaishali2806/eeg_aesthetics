%% Define some stuff!
tic
clc %Remove old stuff
clear all
close all
load('D:\Dropbox\Internship\Stimilu task\data-ort.mat')
try
    
    start_time=GetSecs; %Start timer
    
    rng('shuffle'); %Set random numbers to be truly random (this can throw an error -> remove it then, it does the same more clumsily below)
    %rng(10)
    addpath('D:\Dropbox\Internship\Stimilu task\art_all_renamed\'); %Add the path to stimuli
    
    cfg.eeg_mode=0; %toggle EEG mode (set to 0 to run it without the EEG)
    
    cfg.windowcolor=[125,125,125]; %Colors and Fonts
    cfg.textcolor=[255,255,255];
    cfg.fixcolor=[225,25,225];
    cfg.fontname='Arial';
    cfg.fontsize=20;
    
    cfg.nStim=999;
    cfg.nRep=12;
    cfg.trialamount=100; %cfg.nRep*cfg.nStim;  % n Stimuli, n Repetitions
    cfg.blocksize = 50;
    
    cfg.picsize_hor=300; %Size of the images (this needs to be square for the ratings, pix size is fixed below)
    cfg.picsize_vert=300;
    
    cfg.deg_object=10; %desired degree (visual angle) of the object
    
    [w,h]=Screen('DisplaySize',0); %Get display size in mm
    cfg.screen_hor=w/10; %Convert to cm (horizontal)
    cfg.screen_vert=h/10; %Convert to cm (vertial)
    
    cfg.screen_dist=60; %Distance to screen in cm
    
    w=Screen('Resolution',0); %Get screen properties
    cfg.screensize=[w.width,w.height]; %Store horizontal/vertical pixels (i.e., resolution)
    
    cfg.time_stim=.15; %Duration of the Stimulus
    cfg.time_wait=1; %quick wait between after the stimulus
    cfg.time_iti=.15; %Duration of the inter-trial interval (will be jittered by +/-200ms)
    
    KbName('UnifyKeyNames'); %Fix Mac/Windows differences (hopefully)
    cfg.escapeKey=KbName('q'); %Get key Codes (indices in the KeyCode Matrix) for a few keys
    cfg.spaceKey=KbName('space');
    
    dat.subjcode=input('Enter subject number: '); %Ask for subject number
    
    
    
    
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
    
    real_iti=1000*cfg.time_iti; %For the stimulus
    iti_frames=real_iti/(1000/cfg.monitorFrameRate);
    time_iti=cfg.monitorFlipInterval*(iti_frames-0.5);
    
    clear hidecursor;
    %HideCursor;
    
    screensize=cfg.screensize;
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
    
    %%
    num_blocks = round(cfg.nStim * cfg.nRep / cfg.blocksize);
    
    %num_distractor = num_blocks * 4;
    
    %distractors = randi(6,1,num_blocks);
    
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
    
                                dat = add_distractors(dat, cfg.blocksize);
    
    dat = generate_mask_indices(dat,200);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Start Presentation
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    DrawFormattedText(mainwindow,['instructions \n \n',...
        'Please keep central fixation and don''t blink during the image presentation.\n \n ', ...
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
        img=imread(dat.new_stim{trial});
        img_rand_orig = imread(['Masks_' num2str(dat.mask_indices(trial)) '.jpg']);
        img_rand = desaturate_image1(img_rand_orig);
        stim_tex_rand=Screen('MakeTexture',mainwindow,img_rand);
        stim_tex=Screen('MakeTexture',mainwindow,img);
        Screen('DrawTexture',mainwindow,stim_tex_rand,[],stimrectangle);
        
        if dat.new_ort{trial} == 1
            
            rect_size = [(stimrectangle(3) - stimrectangle(1)) (stimrectangle(4) - stimrectangle(2)) ];
            img_size = rect_size(2) * dat.new_asp{trial};
            diff = rect_size(1) - img_size;
            Screen('DrawTexture',mainwindow,stim_tex,[],[(stimrectangle(1) + diff/2)  stimrectangle(2) (stimrectangle(1)+ diff/2 + img_size)  stimrectangle(4)]);
            
        else
            rect_size = [(stimrectangle(3) - stimrectangle(1)) (stimrectangle(4) - stimrectangle(2)) ]
            img_size = rect_size(1) / dat.new_asp{trial}
            diff = rect_size(2) - img_size;
            Screen('DrawTexture',mainwindow,stim_tex,[],[stimrectangle(1)   (stimrectangle(2)+ diff/2) stimrectangle(3)  (stimrectangle(2) + img_size + diff/2)]);
            
        end
        Screen('DrawDots',mainwindow,fixposition,5,cfg.fixcolor);
         if trial == 1
            ON =Screen('Flip',mainwindow);
        else
            ON =Screen('Flip',mainwindow,ON + time_iti);
        end
        
        flip_1_timing = ON;
        
        if trial ~= 1
            flip_duration(trial,2) = flip_1_timing - flip_2_timing;
        end
        
        
        
        %stimulus off
        Screen('DrawDots',mainwindow,fixposition,5,cfg.fixcolor);
        ON=Screen('Flip',mainwindow,ON + time_stim);
        
        flip_2_timing =ON;
        flip_duration(trial,1) = flip_2_timing - flip_1_timing;
        %
        
        if ~(rem(trial, cfg.blocksize))
            
            a1=360/2;
            start_angle=360*rand;  % A random start
            cm=colormap(gray(4))*255;
            for i=1:2
                a0=start_angle+(i-1)*a1;
                Screen('FillArc',mainwindow,cm(i+2,:),wheelrectangle,a0,a1);
                r=1.25*(wheelrectangle(3)-wheelrectangle(1))/2;
                x_text=fixposition(1)+r*cos(deg2rad(start_angle+(i-1)*a1-1.5*a1-90));
                y_text=fixposition(2)+r*sin(deg2rad(start_angle+(i-1)*a1-1.5*a1-90));
                if i==1
                    DrawFormattedText(mainwindow,'no',x_text,y_text,cfg.textcolor);
                else
                    DrawFormattedText(mainwindow,'yes',x_text,y_text,cfg.textcolor);
                end
            end
            Screen('FillOval',mainwindow,cfg.windowcolor,picrectangle);
            DrawFormattedText(mainwindow,'Found the Pikachu art ?','center','center',cfg.textcolor);
            ON=Screen('Flip',mainwindow,ON+cfg.time_wait);
            
            %show cursor
            SetMouse(fixposition(1),fixposition(2));
            ShowCursor('Hand');
            
            
            %keep waiting for them to select something
            dat.resp1(trial)=0;
            cm=colormap(gray(9))*255;
            while dat.resp1(trial)==0 || keyCode(cfg.spaceKey)~=1
                [keyIsDown,timeSecs,keyCode]=KbCheck(-1);
                [x,y,buttons]=GetMouse;
                if buttons(1)~=0
                    for i=1:2
                        r=max(cfg.screensize);
                        ax0=deg2rad(start_angle+(i-1)*a1-90);
                        ax1=deg2rad(start_angle+i*a1-90);
                        if ax0>ax1
                            t=linspace(ax1,ax0);
                        else
                            t=linspace(ax0,ax1);
                        end
                        x_arc=[fixposition(1),fixposition(1)+r*cos(t),fixposition(1)];
                        y_arc=[fixposition(2),fixposition(2)+r*sin(t),fixposition(2)];
                        cursor_dist=dist([x,y;fixposition(1),fixposition(2)]');
                        if inpolygon(x,y,x_arc,y_arc)==1 && cursor_dist(2)>=(picrectangle(3)-picrectangle(1))/2 && cursor_dist(2)<=(wheelrectangle(3)-wheelrectangle(1))/2
                            for j=1:2
                                a0=start_angle+(j-1)*a1;
                                if j==i
                                    Screen('FillArc',mainwindow,cfg.fixcolor,wheelrectangle,a0,a1);
                                else
                                    Screen('FillArc',mainwindow,cm(j+2,:),wheelrectangle,a0,a1);
                                end
                                r=1.25*(wheelrectangle(3)-wheelrectangle(1))/2;
                                x_text=fixposition(1)+r*cos(deg2rad(start_angle+(j-1)*a1-1.5*a1-90));
                                y_text=fixposition(2)+r*sin(deg2rad(start_angle+(j-1)*a1-1.5*a1-90));
                                if j==1
                                    DrawFormattedText(mainwindow,'no',x_text,y_text,cfg.textcolor);
                                else
                                    DrawFormattedText(mainwindow,'yes',x_text,y_text,cfg.textcolor);
                                end
                            end
                            Screen('FillOval',mainwindow,cfg.windowcolor,picrectangle);
                            DrawFormattedText(mainwindow,'Found the Pikachu art ?','center','center',cfg.textcolor);
                            Screen('Flip', mainwindow);
                            dat.resp1(trial)=i;
                        end
                    end
                end
              
            end
            
            
            %question
            %     a1=360/8;
            %     start_angle=360*rand;  % A random start
            %
            %     for i=1:7
            %         a0=start_angle+(i-1)*a1;
            %         Screen('FillArc',mainwindow,cm(i+2,:),wheelrectangle,a0,a1);
            %         r=1.2*(wheelrectangle(3)-wheelrectangle(1))/2;
            %         x_text=fixposition(1)+r*cos(deg2rad(start_angle+(i-1)*a1-1.5*a1));
            %         y_text=fixposition(2)+r*sin(deg2rad(start_angle+(i-1)*a1-1.5*a1));
            %         DrawFormattedText(mainwindow,num2str(i),x_text,y_text,cfg.textcolor);
            %     end
            %     Screen('FillOval',mainwindow,cfg.windowcolor,picrectangle);
            %     DrawFormattedText(mainwindow,'how many times you found a pikachu in the block?','center','center',cfg.textcolor);
            %     ON=Screen('Flip',mainwindow,ON+ cfg.time_wait);
            %
            %     %keep waiting for them to select something
            %     dat.resp2(trial)=0;
            %     while dat.resp2(trial)==0 || keyCode(cfg.spaceKey)~=1
            %         [keyIsDown,timeSecs,keyCode]=KbCheck(-1);
            %         [x,y,buttons]=GetMouse;
            %         if buttons(1)~=0
            %             for i=1:7
            %                 r=max(cfg.screensize);
            %                 ax0=deg2rad(start_angle+(i-1)*a1-90);
            %                 ax1=deg2rad(start_angle+i*a1-90);
            %                 if ax0>ax1
            %                     t=linspace(ax1,ax0);
            %                 else
            %                     t=linspace(ax0,ax1);
            %                 end
            %                 x_arc=[fixposition(1),fixposition(1)+r*cos(t),fixposition(1)];
            %                 y_arc=[fixposition(2),fixposition(2)+r*sin(t),fixposition(2)];
            %                 cursor_dist=dist([x,y;fixposition(1),fixposition(2)]');
            %                 if inpolygon(x,y,x_arc,y_arc)==1 && cursor_dist(2)>=(picrectangle(3)-picrectangle(1))/2 && cursor_dist(2)<=(wheelrectangle(3)-wheelrectangle(1))/2
            %                     for j=1:7
            %                         a0=start_angle+(j-1)*a1;
            %                         if j==i
            %                             Screen('FillArc',mainwindow,cfg.fixcolor,wheelrectangle,a0,a1);
            %                         else
            %                             Screen('FillArc',mainwindow,cm(j+2,:),wheelrectangle,a0,a1);
            %                         end
            %                         r=1.2*(wheelrectangle(3)-wheelrectangle(1))/2;
            %                         x_text=fixposition(1)+r*cos(deg2rad(start_angle+(j-1)*a1-1.5*a1));
            %                         y_text=fixposition(2)+r*sin(deg2rad(start_angle+(j-1)*a1-1.5*a1));
            %                         DrawFormattedText(mainwindow,num2str(j),x_text,y_text,cfg.textcolor);
            %                     end
            %                     Screen('FillOval',mainwindow,cfg.windowcolor,picrectangle);
            %                     DrawFormattedText(mainwindow,'how many times you found a pikachu in the block?','center','center',cfg.textcolor);
            %                     Screen('Flip', mainwindow);
            %                     dat.resp2(trial)=i;
            %                 end
            %             end
            %        end
            %    end
            
        end
        
        HideCursor;
        
        %KbWait(-1);
        
        %close texture
        Screen('Close',stim_tex_rand);
        
        
        
        
    end
    
    sca
    save(['RSVP_eeg_s',num2str(dat.subjcode)],'cfg','dat');
    
catch
    
   sca
   lasterr
    
end

toc


