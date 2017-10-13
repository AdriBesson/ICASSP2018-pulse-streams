%% Monte Carlo simulation - 1 channel
clear all;
close all;
clc
addpath(genpath('src'));
addpath(genpath('src/utils'));

%-- Generate the US sequence object
us_seq = USSequence();

%-- Number of pulses in each channel
n_pulses = 20;

%-- Noise level
noise_level = 1000;

%-- Number of simulation draws
n_draws = 1000;

%-- Measurement ratios
meas_ratio = [0.01:0.01:0.20, 0.30, 0.40, 0.50];

%-- Number of channels as prior
n_cha_prior = 0;

%-- L1 regularization - inter-element spacing 1 wavelength
[nmse, nrmse] = benchmark(us_seq, 'l1', meas_ratio, n_cha_prior, n_pulses, n_draws, noise_level);

%-- Save the output file
filenameOut = 'results/results_1channel_synth_pulse.mat';
save(filenameOut, 'nmse', 'nrmse', 'us_seq', 'n_pulses', 'n_draws', 'meas_ratio');

%-- Least squares - inter-element spacing 1 wavelength
[nmse, nrmse] = benchmark(us_seq, 'LS', meas_ratio, n_cha_prior, n_pulses, n_draws, noise_level);

%-- Save the output file
filenameOut = 'results/results_1channel_least_squares_synth_pulse.mat';
save(filenameOut, 'nmse', 'nrmse', 'us_seq', 'n_pulses', 'n_draws', 'meas_ratio');

%-- L1 regularization - inter-element spacing 2 wavelengths
us_seq.set_element_spacing(2*us_seq.speed_of_sound/us_seq.central_frequency);
[nmse, nrmse] = benchmark(us_seq, 'l1', meas_ratio, n_cha_prior, n_pulses, n_draws, noise_level);

%-- Save the output file
filenameOut = 'results/results_1channel_synth_pulse_2lambda.mat';
save(filenameOut, 'nmse', 'nrmse', 'us_seq', 'n_pulses', 'n_draws', 'meas_ratio');

%-- Least squares - inter-element spacing 2 wavelengths
us_seq.set_element_spacing(2*us_seq.speed_of_sound/us_seq.central_frequency);
[nmse, nrmse] = benchmark(us_seq, 'LS', meas_ratio, n_cha_prior, n_pulses, n_draws, noise_level);

%-- Save the output file
filenameOut = 'results/results_1channel_least_squares_synth_pulse_2lambda.mat';
save(filenameOut, 'nmse', 'nrmse', 'us_seq', 'n_pulses', 'n_draws', 'meas_ratio');

%% Monte Carlo simulation - multiple channels
clear all;
close all;
clc
addpath(genpath('../src'));
addpath(genpath('../src/utils'));
addpath(genpath('../data'));

%-- Generate the US sequence object
us_seq = USSequence();
xm = us_seq.get_element_locations();

%-- Number of pulses in each channel
n_pulses = 20;

%-- Noise level
noise_level = 1000; % no noise

%-- Number of simulation draws
n_draws = 1000;

%-- Measurement ratios
meas_ratio = [0.01:0.01:0.20, 0.30, 0.40, 0.50];

%-- List of number of channels onto which prior knowledge is known
list_n_cha_prior = [1, 4, 9, 19];

%-- L1 regularization - inter-element spacing 1 wavelength 
for jj = 1:numel(list_n_cha_prior)
    n_cha_prior = list_n_cha_prior(jj);
    
    % Benchmark
    [nmse, nrmse] = benchmark(us_seq, 'l1', meas_ratio, n_cha_prior, n_pulses, n_draws, noise_level);
    
    %-- Save the output file
    filenameOut = strcat(['results/','results_', num2str(n_cha_prior+1) ,'channels_synth_pulse.mat']);
    save(filenameOut, 'nmse', 'nrmse', 'us_seq', 'n_pulses', 'n_draws', 'meas_ratio');
end

%-- Least-squares - inter-element spacing 1 wavelength
for jj = 1:numel(list_n_cha_prior)
    n_cha_prior = list_n_cha_prior(jj);
    
    % Benchmark
    [nmse, nrmse] = benchmark(us_seq, 'LS', meas_ratio, n_cha_prior, n_pulses, n_draws, noise_level);
    
    %-- Save the output file
    filenameOut = strcat(['results/','results_', num2str(n_cha_prior+1) ,'channels_least_squares_synth_pulse.mat']);
    save(filenameOut, 'nmse', 'nrmse', 'us_seq', 'n_pulses', 'n_draws', 'meas_ratio');
end

%-- L1 regularization - inter-element spacing 2 wavelengths
us_seq.set_element_spacing(2*us_seq.speed_of_sound/us_seq.central_frequency);
for jj = 1:numel(list_n_cha_prior)
    n_cha_prior = list_n_cha_prior(jj);
    
    % Benchmark
    [nmse, nrmse] = benchmark(us_seq, 'l1', meas_ratio, n_cha_prior, n_pulses, n_draws, noise_level);
    
    %-- Save the output file
    filenameOut = strcat(['results/','results_', num2str(n_cha_prior+1) ,'channels_synth_pulse_2lambda.mat']);
    save(filenameOut, 'nmse', 'nrmse', 'us_seq', 'n_pulses', 'n_draws', 'meas_ratio');
end

%-- Least-squares - inter-element spacing 1 wavelength
us_seq.set_element_spacing(2*us_seq.speed_of_sound/us_seq.central_frequency);
for jj = 1:numel(list_n_cha_prior)
    n_cha_prior = list_n_cha_prior(jj);
    
    % Benchmark
    [nmse, nrmse] = benchmark(us_seq, 'l1', meas_ratio, n_cha_prior, n_pulses, n_draws, noise_level);
    
    %-- Save the output file
    filenameOut = strcat(['results/','results_', num2str(n_cha_prior+1) ,'channels_least_squares_synth_pulse_2lambda.mat']);
    save(filenameOut, 'nmse', 'nrmse', 'us_seq', 'n_pulses', 'n_draws', 'meas_ratio');
end