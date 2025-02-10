
function data = generate_mask_indices(dat, unique_ent)

length_vec = length(dat.new_stim)
% Set the length of the vector and number of unique entries
vector_length = length_vec;
unique_entries = unique_ent;

% Generate a vector with equal entries for each number from 1 to 20
random_vector = repmat(1:unique_entries, 1, round(vector_length / unique_entries));

% Shuffle the vector to randomize the order
random_vector = random_vector(randperm(length(random_vector)));

% Display the generated vector
disp(random_vector);
dat.mask_indices = random_vector;
data = dat

end 