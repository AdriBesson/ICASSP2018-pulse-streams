function [nmse, nrmse] = benchmark(us_seq, meas_ratio, n_cha_prior, n_points, n_draws)
nmse = zeros(numel(meas_ratio), n_draws);
nrmse = zeros(numel(meas_ratio), n_draws);

%-- Setting the sparsity model
sparsity_model = SparsityModel();
pulse = us_seq.estimate_received_pulse();

%-- Simulation
for mm = 1:numel(meas_ratio)
    M = round(meas_ratio(mm)*us_seq.number_time_samples);
    for kk = 1:n_draws
        
        % positions of the points
        points_locations(:,1) = randi(us_seq.number_elements, [n_points, 1]);
        points_locations(:,2) = randi([floor(numel(pulse)/2) us_seq.number_time_samples-floor(numel(pulse)/2)], [n_points, 1]);
        points_amplitudes =  rand(size(points_locations(:,2)));
        
        %-- Positions and amplitude of the points
        [raw_data, points_locations_raw] = us_seq.generate_rawdata(points_locations, points_amplitudes, 100);
        
        % Reference channels
        channel_number_prior = randi(us_seq.number_elements-n_cha_prior-1) + (1:n_cha_prior+1);
        
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
                for ll = 1:n_points
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
        
        % solving the problem
        alpha_est = solver.solve();
        channel_est = us_seq.get_pulse_dictionary()*alpha_est;
        
        % computing the nmse and nrmse
        nmse(mm,kk) = norm(channel_est - channel)^2 ./ norm(channel)^2;
        nrmse(mm,kk) = goodnessOfFit(channel_est, channel, 'NRMSE');
        
        % log
        fprintf('************** Experiment %i, Measurement ratio %e **************\n',kk, meas_ratio(mm));
        fprintf(' NMSE = %e, NRMSE = %e\n', ...
            nmse(mm,kk), nrmse(mm,kk));
    end
end

