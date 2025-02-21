
function data_out = add_distractors_faster(dat, blocksize)
%%% This function adds triggers

dat = dat;
blocksize = blocksize; 
num_dist = length(dat.stim)/blocksize; 

% generate a sequence of randomised 0 & 1 

dist_block = randi([0, 1], [1, round(num_dist)]); 

n = round(num_dist); % 
numberOfOnes = n/2;
indexes = randperm(n);
dist_block = zeros(1, n);
dist_block(indexes(1:numberOfOnes)) = 1;

% create vector having n distractor images repated 6 times 
dist_ind = []
     for nn=1:24
    dist_ind=[dist_ind;randperm(5)'];
end

dist_img_blk = 0;
dat.new_stim = {};
dat.new_asp = {};
dat.new_ort = {};


addpath('C:\Masters Cyber\Master Thesis\Thesis_VG\data-ort.mat');

for ii = 1 : length(dist_block); 
    
    if ii < length(dist_block)
         data = dat.stim(1,(ii -1)* 200 + 1 : (ii)* 200);
         aspect_ratio = dat.aspect_ratio(1,(ii -1)* 200 + 1 : (ii)* 200);
         ort = dat.ort(1,(ii -1)* 200 + 1 : (ii)* 200);
    else 
        data = dat.stim(1,(ii -1)* 200 + 1 : end);
        aspect_ratio = dat.ort(1,(ii -1)* 200 + 1 : end);
        ort = dat.ort(1,(ii -1)* 200 + 1 : end);
    end 
    
    
    if dist_block(ii) == 1
        dist_img_blk = dist_img_blk + 1;
        A = data;
        b = ['dist_' num2str(dist_ind(dist_img_blk)) '.jpg'];
        k = randi([1 199]); %row position, can be 0,1,2 or 3 in this case
        data = [A(1,1:k) b A(1,k+1:end)];
        
        % load image
        img = imread(b);
        aspect_rat = size(img,2)/ size(img,1);
        if aspect_rat > 1
            ort_n = 2; 
        else 
            ort_n = 1; 
        end 
    
        A_ort = ort;
        b_ort = ort_n;%['dist_' num2str(dist_ind(dist_img_blk)) '.jpg'];
        data_ort = [A_ort(1,1:k) b_ort A_ort(1,k+1:end)];
        
        
        A_asp = aspect_ratio;
        b_asp = aspect_rat;%['dist_' num2str(dist_ind(dist_img_blk)) '.jpg'];
        data_asp = [A_asp(1,1:k) b_asp A_asp(1,k+1:end)];
        
        
        
        
    else 
        data = data;
        data_ort = ort;
        data_asp = aspect_ratio;
    end 
    
    dat.new_stim = [dat.new_stim data];
    dat.new_asp = [dat.new_asp data_asp];
    dat.new_ort = [dat.new_ort data_ort];
    dat.dist_block = dist_block;
    
    data_out = dat;
    
end 

end 
