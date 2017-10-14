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

%-- Generate the US sequence object
us_seq = USSequence();
xm = us_seq.get_element_locations();

%-- Number of pulses in each channel
n_pulses = 20;

%-- Noise level
noise_level = 30; % noise level at -40dB
epsilon = sqrt(10^(-(noise_level)/10)); % Bound depending on the noise level
regularization_parameter = 5e0;

%-- Measurement ratios
meas_ratio = 0.08;

%-- Generate the raw data
pulse = us_seq.estimate_received_pulse();

% positions of the points
points_locations(:,1) = randi(us_seq.number_elements, [n_pulses, 1]);
points_locations(:,2) = randi([floor(numel(pulse)/2) us_seq.number_time_samples-floor(numel(pulse)/2)], [n_pulses, 1]);
points_amplitudes =  rand(size(points_locations(:,2)));
points_locations_raw = us_seq.generate_rawdata(points_locations, points_amplitudes, 1000);
raw_data = us_seq.data;

%-- Reference channel
channel_number = 32;
channel = raw_data(:,channel_number);

%-- Noisy raw data
raw_data_noisy = awgn(raw_data, noise_level);
channel_noisy = raw_data_noisy(:,channel_number);

%% Single channel experiment
n_cha_prior = 0; 

%-- Reconstruction
channel_est = reconstruct_image(us_seq, meas_ratio, raw_data_noisy, channel_number, n_cha_prior, points_locations_raw, epsilon, regularization_parameter);

%-- Save the output file
t = us_seq.get_time_samples();
filenameOut = 'results/results_1channel_noisy.mat';
save(filenameOut, 't', 'channel_noisy', 'channel_est', 'channel');

%% Multi-channel experiment
n_cha_prior = 4; 

%-- Reconstruction
channel_est = reconstruct_image(us_seq, meas_ratio, raw_data_noisy, channel_number, n_cha_prior, points_locations_raw, epsilon, regularization_parameter);

%-- Save the output file
filenameOut = 'results/results_multichannel_noisy.mat';
save(filenameOut, 't', 'channel_noisy', 'channel_est', 'channel');
