%-- Script used to generate Monte Carlo simulation - multi-channel
clear all;
close all;
clc
addpath(genpath('../src'));

%-- Pulse
f0 = 5e6;
fs = 4*1.5*f0;
c0 = 1540;
N_el = 64;
xm = (0:N_el-1)*2*c0/f0;
xm = xm - xm(end/2);
impulse = sin(2*pi*f0*(0:1/fs:2/f0));
impulse_response = impulse .* hanning(length(impulse))';
excitation = sin(2*pi*f0*(0:1/fs:1/f0));
pulse = conv(conv(excitation, impulse_response), impulse_response);
pulse = pulse /max(abs(pulse));

%-- Characteristics of the raw data
Nt = 2000;
n_points = 20;

%-- Build the dictionary
pulse_pad = zeros(Nt,1);
pulse_pad(1:numel(pulse),:) = pulse;
Psi = circulant(pulse_pad);

%-- Number of simulation draws
n_draws = 1000;

%-- Measurement ratios
meas_ratio = [0.01:0.01:0.2, 0.3:0.1:0.5];

%-- Output variables
nmse = zeros(numel(meas_ratio), n_draws);
nrmse = zeros(numel(meas_ratio), n_draws);

%-- List of number of channels onto which prior knowledge is known
list_n_cha_prior = [1, 4, 9, 19];

%-- Number of channels used for the prior
for jj = 1:numel(list_n_cha_prior)
    n_cha_prior = list_n_cha_prior(jj);

    %-- Simulation
    for mm = 1:numel(meas_ratio)
        M = round(meas_ratio(mm)*Nt);
        for kk = 1:n_draws

            % positions of the points
            x_points = randi(N_el, [n_points, 1]);
            t_points = randi([floor(numel(pulse)/2) Nt-floor(numel(pulse)/2)], [n_points, 1]);
            amp_points =  rand(size(t_points));

            % Creation of the raw data - 1PW insonification
            point_locations = zeros(N_el, n_points);
            rf_data = zeros(Nt, N_el);
            for pp = 1:n_points
                for ll = 1:N_el
                    t_Tx = t_points(pp) / fs / 2;
                    t_Rx = sqrt(1/c0^2*(xm(x_points(pp))-xm(ll))^2 + (t_points(pp) / 2 / fs)^2);
                    pos_raw_points = max(round((t_Tx + t_Rx)*fs-1), 1);
                    pos_raw_points = min(pos_raw_points, Nt);
                    point_locations(ll,pp) = pos_raw_points;
                    rf_data( pos_raw_points ,ll) = amp_points(pp);
                end
            end
            for ll = 1:N_el
                rf_data_conv(:,ll) = conv(rf_data(:,ll), pulse);
            end
            rf_data = rf_data_conv(1:Nt,:);

            % Reference channels
            channel_number_prior = randi(N_el-n_cha_prior-1) + (1:n_cha_prior+1);

            % Considered channel
            channel_number = channel_number_prior(round(numel(channel_number_prior)/2));
            channel = rf_data(:,channel_number);
            channel = channel ./ max(channel);
            channel_number_prior = channel_number_prior(not(channel_number_prior == channel_number));

            % Creation of the mask for the support
            mask = ones(numel(channel),1);
            for pp = 1:numel(channel_number_prior)
                mask_pp = zeros(numel(channel),1);
                for ll = 1:n_points
                    delta_ind = round(abs(xm(channel_number_prior(pp)) - xm(channel_number))*fs / c0);
                    cur_point = point_locations(channel_number_prior(pp), ll);
                    mask_tmp = zeros(size(mask_pp));
                    mask_tmp(max(cur_point-delta_ind, 1):min(cur_point+delta_ind, numel(mask_pp))) = 1;
                    mask_pp = mask_pp + mask_tmp;
                end
                mask = and(mask, mask_pp);
            end
            mask(mask > 1) = 1;

            %-- Measurement matrix and measurement vector
            Phi = randraw('normal', [0, 1/sqrt(M)], [M Nt]);
            y = Phi*channel(:);

            %-- Setting the functions for the admm solver
            G = Phi*Psi;
            A =@(x) G*(mask.*x(:));
            At =@(x) mask.*(G'*x(:));
            T = @(x) x;
            Tt = @(x) x;

            %-- Setting the parameters for the solver
            param_solver.verbose = 0; % display parameter
            param_solver.max_iter = 1000;       % maximum iteration
            param_solver.tol = 1e-10;        % tolerance to stop iterating
            param_solver.nu = pow_method(A,At,[numel(channel), 1],1e-4,100,0);
            param_solver.gamma = 8e-1 *norm(At(y), inf);
            epsilon = 0;

            % solving the problem
            alpha_est = admm_bpcon(y, epsilon, A, At, T, Tt, param_solver);
            channel_est = Psi*alpha_est;

            % computing the nmse and nrmse
            nmse(mm,kk) = norm(channel_est - channel)^2 ./ norm(channel)^2;
            nrmse(mm,kk) = goodnessOfFit(channel_est, channel, 'NRMSE');

            % log
            fprintf('************** Experiment %i, Measurement ratio %e **************\n',kk, meas_ratio(mm));
            fprintf(' NMSE = %e, NRMSE = %e\n', ...
                nmse(mm,kk), nrmse(mm,kk));
        end
    end

    %-- Save the output file
    filenameOut = strcat(['../resultsSPL/','results_', num2str(n_cha_prior+1) ,'channels_synth_pulse_2lambda.mat']);
    save(filenameOut, 'nmse', 'nrmse', 'Nt', 'n_points', 'n_draws', 'f0', 'fs', 'c0', 'pulse', 'xm', 'meas_ratio');
end