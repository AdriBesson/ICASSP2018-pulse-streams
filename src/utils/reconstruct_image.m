function channel_est = reconstruct_image(us_seq, meas_ratio, raw_data, channel_number, n_cha_prior, points_locations_raw, noise_level)

%-- Setting the sparsity model
sparsity_model = SparsityModel();

%-- Simulation
M = round(meas_ratio*us_seq.number_time_samples);

% Reference channels
channel_number_prior = channel_number + (-n_cha_prior/2:n_cha_prior/2);

% Considered channel
channel_number = channel_number_prior(round(numel(channel_number_prior)/2));
channel = raw_data(:,channel_number);
channel = channel ./ max(channel);
channel_number_prior = channel_number_prior(not(channel_number_prior == channel_number));

%-- Measurement matrix and measurement vector
Phi = randraw('normal', [0, 1/sqrt(M)], [M, us_seq.number_time_samples]);
y = Phi*channel(:);

%-- Setting the functions for the admm solver
measurement_matrix = Phi*us_seq.get_pulse_dictionary();
measurement_model = MeasurementModel(measurement_matrix);

% Creation of the mask for the support
mask = ones(numel(channel),1);
if (~isempty(channel_number_prior))
    xm = us_seq.get_element_locations();
    for pp = 1:numel(channel_number_prior)
        mask_pp = zeros(numel(channel),1);
        for ll = 1:size(points_locations_raw, 2)
            delta_ind = round(abs(xm(channel_number_prior(pp)) - xm(channel_number))*us_seq.sampling_frequency / us_seq.speed_of_sound);
            cur_point = points_locations_raw(channel_number_prior(pp), ll);
            mask_tmp = zeros(size(mask_pp));
            mask_tmp(max(cur_point-delta_ind, 1):min(cur_point+delta_ind, numel(mask_pp))) = 1;
            mask_pp = mask_pp + mask_tmp;
        end
        mask = and(mask, mask_pp);
    end
    mask(mask > 1) = 1;
    measurement_model.set_mask(mask);
end

%-- Setting the parameters for the solver
solver = Solver(measurement_model, sparsity_model, y);

%-- Change the radius of the l2-ball depending on the noise level
solver.set_radius(sqrt(10^(-(noise_level)/10)));

%-- Change the regularization parameter
solver.set_regularization_parameter(5e0 *norm(measurement_model.backward(y), inf));

% solving the problem
alpha_est = solver.solve();
channel_est = us_seq.get_pulse_dictionary()*alpha_est;

end

