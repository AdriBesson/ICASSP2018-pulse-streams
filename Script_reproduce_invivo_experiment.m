%-- Script used to generate Monte Carlo simulation - 1 channel
clear all;
close all;
clc

addpath(genpath('src'));
addpath(genpath('src/utils'));
addpath(genpath('data'));

%-- Fix the seed number for reproducibility
rng('default');
rng(1);

%-- Load the data file
load('data/invivo_CND_data.mat')

%-- Generate the US sequence object
us_seq = USSequence();
us_seq_rec_1channel = USSequence();
us_seq_rec_2channels = USSequence();
rawdata = rawdata(50:900,:);

%-- Set the parameters of the sequence 
us_seq.set_bandwidth(1);
us_seq.set_sampling_frequency(setup.Fs*1e6);
us_seq.set_central_frequency(setup.Fcenter*1e6);
us_seq.set_number_elements(setup.Nex_p1);
us_seq.set_number_time_samples(size(rawdata, 1));
us_seq.set_element_spacing(setup.dx_p1*1e-3);
us_seq.set_speed_of_sound(setup.c*1e3);
us_seq.set_excitation_cycle(0.5);
us_seq.set_initial_time(50/us_seq.sampling_frequency);
us_seq.set_impulse_response_cycle(1);
us_seq.set_data(rawdata);

us_seq_rec_1channel.set_bandwidth(1);
us_seq_rec_1channel.set_sampling_frequency(setup.Fs*1e6);
us_seq_rec_1channel.set_central_frequency(setup.Fcenter*1e6);
us_seq_rec_1channel.set_number_elements(setup.Nex_p1);
us_seq_rec_1channel.set_number_time_samples(size(rawdata, 1));
us_seq_rec_1channel.set_element_spacing(setup.dx_p1*1e-3);
us_seq_rec_1channel.set_speed_of_sound(setup.c*1e3);
us_seq_rec_1channel.set_excitation_cycle(0.5);
us_seq_rec_1channel.set_impulse_response_cycle(1);
us_seq_rec_1channel.set_data(rawdata);
us_seq_rec_1channel.set_initial_time(50/us_seq.sampling_frequency);

us_seq_rec_2channels.set_bandwidth(1);
us_seq_rec_2channels.set_sampling_frequency(setup.Fs*1e6);
us_seq_rec_2channels.set_central_frequency(setup.Fcenter*1e6);
us_seq_rec_2channels.set_number_elements(setup.Nex_p1);
us_seq_rec_2channels.set_number_time_samples(size(rawdata, 1));
us_seq_rec_2channels.set_element_spacing(setup.dx_p1*1e-3);
us_seq_rec_2channels.set_speed_of_sound(setup.c*1e3);
us_seq_rec_2channels.set_excitation_cycle(0.5);
us_seq_rec_2channels.set_impulse_response_cycle(1);
us_seq_rec_2channels.set_data(rawdata);
us_seq_rec_2channels.set_initial_time(50/us_seq.sampling_frequency);

%-- Measurement ratio
meas_ratio = 0.03;


%% 1-channel scenario

%-- Number of prior channels
n_cha_prior = 0;

%-- Reconstruction
regularization_parameter =  3e-1;
rawdata_rec_1channel = zeros(size(rawdata));
disp('************ In-vivo data - 1-channel scenario************')
for kk = 1:us_seq.number_elements
    rawdata_rec_1channel(:,kk) = reconstruct_image(us_seq_rec_1channel, meas_ratio, rawdata, kk, n_cha_prior, [], [], regularization_parameter);
end
us_seq_rec_1channel.set_data(rawdata_rec_1channel);

%% 2-channel scenario
rawdata_rec_2channels = zeros(size(rawdata));

%-- Initialization: reconstruction of the reference channel
n_cha_prior = 0;
channel_number = 15;
[rawdata_rec_2channels(:,channel_number), alpha_est_ref] = reconstruct_image(us_seq_rec_2channels, 1, rawdata, channel_number, n_cha_prior, [], [], regularization_parameter);

%-- Number of prior channels
n_cha_prior = 1;

%-- Reconstruction
number_of_points_support = 30;
regularization_parameter =  3e-1;
alpha_est = alpha_est_ref;
disp('************ In-vivo data - 2-channel scenario************')
for kk = channel_number-1:-1:1
    points_locations_raw = get_support(us_seq, alpha_est, number_of_points_support);
    [rawdata_rec_2channels(:,kk), alpha_est] = reconstruct_image(us_seq_rec_2channels, meas_ratio, rawdata, kk, n_cha_prior, points_locations_raw', [], regularization_parameter);
end
alpha_est = alpha_est_ref;
for kk = channel_number+1:us_seq.number_elements
    points_locations_raw = get_support(us_seq, alpha_est, number_of_points_support);
    [rawdata_rec_2channels(:,kk), alpha_est] = reconstruct_image(us_seq_rec_2channels, meas_ratio, rawdata, kk, n_cha_prior, points_locations_raw', [], regularization_parameter);
end
us_seq_rec_2channels.set_data(rawdata_rec_2channels);

%% Beamforming, postprocessing and saving
%-- Image grid
x = us_seq.get_element_locations();
t = us_seq.get_time_samples();
z = us_seq.speed_of_sound*t / 2;

xim = x(1):us_seq.speed_of_sound/us_seq.central_frequency/4:10/1000;
zim = 18/1000:us_seq.speed_of_sound/us_seq.central_frequency/8:z(end);

%-- Beamforming and post-processing
bmode_ref = postprocess(us_seq.beamform(xim, zim));
bmode_rec_1channel = postprocess(us_seq_rec_1channel.beamform(xim, zim));
bmode_rec_2channels = postprocess(us_seq_rec_2channels.beamform(xim, zim));

%-- Saving
filenameOut = 'resultsSPL/results_invivo_cnd.mat';
save(filenameOut, 'bmode_ref', 'bmode_rec_1channel', 'bmode_rec_2channels', 'xim', 'zim')

