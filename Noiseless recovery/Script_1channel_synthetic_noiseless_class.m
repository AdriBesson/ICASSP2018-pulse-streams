%-- Script used to generate Monte Carlo simulation - 1 channel
clear all;
close all;
clc
addpath(genpath('../src'));
addpath(genpath('../src/utils'));
addpath(genpath('../data'));

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

%% inter-element spacing 1 wavelength
%-- Simulation
[nmse, nrmse] = benchmark(us_seq, meas_ratio, n_cha_prior, n_pulses, n_draws, noise_level);

%-- Save the output file
filenameOut = '../resultsSPL/results_1channel_synth_pulse.mat';
save(filenameOut, 'nmse', 'nrmse', 'us_seq', 'n_pulses', 'n_draws', 'meas_ratio');

%% inter-element spacing 2 wavelengths
%-- Simulation
us_seq.set_element_spacing = 2*us_seq.speed_of_sound/us_seq.central_frequency;
[nmse, nrmse] = benchmark(us_seq, meas_ratio, n_cha_prior, n_pulses, n_draws, noise_level);

%-- Save the output file
filenameOut = '../resultsSPL/results_1channel_synth_pulse_2lambda.mat';
save(filenameOut, 'nmse', 'nrmse', 'us_seq', 'n_pulses', 'n_draws', 'meas_ratio');