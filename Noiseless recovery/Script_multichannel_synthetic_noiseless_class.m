%-- Script used to generate Monte Carlo simulation - multi-channel
clear all;
close all;
clc
%-- Script used to generate Monte Carlo simulation - 1 channel
clear all;
close all;
clc
addpath(genpath('src'));
addpath(genpath('src/utils'));
addpath(genpath('data'));

%-- Generate the US sequence object
us_seq = USSequence();
xm = us_seq.get_element_locations();

%-- Number of pulses in each channel
n_pulses = 20;

%-- Number of simulation draws
n_draws = 1000;

%-- Measurement ratios
meas_ratio = [0.01:0.01:0.20, 0.30, 0.40, 0.50];

%-- List of number of channels onto which prior knowledge is known
list_n_cha_prior = [1, 4, 9, 19];

%-- Setting the sparsity model
sparsity_model = SparsityModel();

%-- Number of channels used for the prior
for jj = 1:numel(list_n_cha_prior)
    n_cha_prior = list_n_cha_prior(jj);
    
    % Benchmark
    [nmse, nrmse] = benchmark(us_seq, meas_ratio, n_cha_prior, n_pulses, n_draws);
    
    %-- Save the output file
    filenameOut = strcat(['../resultsSPL/','results_', num2str(n_cha_prior+1) ,'channels_synth_pulse_2lambda.mat']);
    save(filenameOut, 'nmse', 'nrmse', 'Nt', 'n_points', 'n_draws', 'f0', 'fs', 'c0', 'pulse', 'xm', 'meas_ratio');
end