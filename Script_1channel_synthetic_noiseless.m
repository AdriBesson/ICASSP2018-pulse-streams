%-- Script used to generate Monte Carlo simulation - 1 channel
clear all;
close all;
clc
addpath(genpath('src'));
addpath(genpath('data'));
addpath(genpath('../toolboxes/unlocbox'));

%-- Pulse
f0 = 5e6;
fs = 4*1.5*f0;
c0 = 1540;
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
meas_ratio = [0.01:0.01:0.20, 0.30, 0.40, 0.50];

%-- Output variables
nmse = zeros(numel(meas_ratio), n_draws);
nrmse = zeros(numel(meas_ratio), n_draws);

%-- Simulation
for mm = 1:numel(meas_ratio)
    M = round(meas_ratio(mm)*Nt);
    for kk = 1:n_draws
        %-- Positions and amplitude of the points
        t_points = randi(Nt, [n_points, 1]);
        amp_points =  rand(size(t_points));
        
        % Channel
        channel = zeros(Nt,1);
        channel(t_points) = amp_points;
        channel = conv(channel, pulse, 'same');
        
        %-- Measurement matrix and measurement vector
        Phi = randraw('normal', [0, 1/sqrt(M)], [M Nt]);
        y = Phi*channel(:);
        
        %-- Setting the functions for the admm solver
        G = Phi*Psi;
        A =@(x) G*x(:);
        At =@(x) G'*x(:);
        T = @(x) x;
        Tt = @(x) x;
        
        %-- Setting the parameters for the solver
        param_solver.verbose = 0; % display parameter
        param_solver.max_iter = 1000;       % maximum iteration
        param_solver.tol = 1e-10;        % tolerance to stop iterating
        param_solver.nu = eigs(G'*G, 1);
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
filenameOut = 'results_1channel_synth_pulse.mat';
save(filenameOut, 'nmse', 'nrmse', 'Nt', 'n_points', 'n_draws', 'f0', 'fs', 'c0', 'pulse', 'meas_ratio');