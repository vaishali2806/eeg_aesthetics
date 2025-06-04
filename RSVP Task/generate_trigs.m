
function dat_out = generate_trigs(dat)

dat = dat;
trigs = [];

for ii = 1: length(dat.new_stim);

    if contains(dat.new_stim{ii}, 'dist_');
        trigs(ii,1) = 2;
    else
        trigs(ii,1) = 1;
    end
    dat.trigs = trigs; 

    dat_out = dat;
end
end