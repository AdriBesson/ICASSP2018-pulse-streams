%-- Script used to generate Monte Carlo simulation - 1 channel
clear all;
close all;
clc
addpath(genpath('src'));
addpath(genpath('src/utils'));
addpath(genpath('data'));

%-- Generate the US sequence object
us_seq = USSequence();

%-- Number of pulses in each channel
n_pulses = 20;

%-- Number of simulation draws
n_draws = 1000;

%-- Measurement ratios
meas_ratio = [0.01:0.01:0.20, 0.30, 0.40, 0.50];

%-- Number of channels as prior
n_cha_prior = 0;

%-- Simulation
[nmse, nrmse] = benchmark(us_seq, meas_ratio, n_cha_prior, n_pulses, n_draws);

%-- Save the output file
filenameOut = '../resultsSPL/results_1channel_synth_pulse.mat';
save(filenameOut, 'nmse', 'nrmse', 'Nt', 'n_points', 'n_draws', 'f0', 'fs', 'c0', 'pulse', 'meas_ratio');