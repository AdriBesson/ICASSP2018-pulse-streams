%-- Script used to reconstruct a real carotid
clear all;
close all;
clc
addpath(genpath('src'));
addpath(genpath('data'));
addpath(genpath('../toolboxes/export_fig'));

%-- Load the file
load('data/rawdata_carotid_ultrasonix5.mat');
rf_data = rawdata(300:end,:,round(size(rawdata, 3)/2));
nbAngles = numel(info_angle);

probe.dBRange = 40;
probe.fs = 40e6;
probe.f0 = 5e6;
probe.bandwidth = [2.5e6,7.5e6];
probe.c0 = 1540;
probe.Pitch = 0.460/1000;
probe.N_active = size(rawdata,2);
probe.transducer_width = (probe.N_active-1)*probe.Pitch;
probe.xm = linspace(-probe.transducer_width/2,probe.transducer_width/2,probe.N_active)-probe.Pitch/2;  
probe.nbAngles = nbAngles;
probe.angles = 0;

impulse = sin(2*pi*probe.f0*(0:1/probe.fs:2/probe.f0));
impulse_response = impulse .* hanning(length(impulse))';
excitation = sin(2*pi*probe.f0*(0:1/probe.fs:1/probe.f0));
pulse = conv(conv(excitation, impulse_response), impulse_response);
t = (0:size(rf_data, 1)-1)/probe.fs;
probe.N_element = probe.N_active;
probe.bandwidth = [probe.f0*0.5, probe.f0*1.5];
probe.transducerwidth = probe.Pitch*(probe.N_active-1);
probe.xm = linspace(-probe.transducerwidth/2,probe.transducerwidth/2,probe.N_active)-probe.Pitch/2;

%-- Multichannel reconstruction
%-- Build the dictionary
pulse_pad = zeros(size(rf_data, 1), 1);
pulse_pad(1:numel(pulse)) = pulse;
Psi = circulant(pulse_pad);

%% Reconstruction
% Measurement matrix for the multichannel
N = size(rf_data, 1);
M = round(0.06*N);
Phi = randraw('normal', [0, 1/sqrt(M)], [M N]);
G = Phi*Psi;
A =@(x) G*x(:);
At =@(x) G'*x(:);
param_solver.nu = eigs(G'*G,1);
list_Nelements = 1:128;

%% 1 channel 
%-- Reconstructed raw data
rf_data_rec1 = zeros(size(rf_data));

for kk = 1:numel(list_Nelements)
    disp(['Reconstruction of channel ', num2str(list_Nelements(kk))]);
    
    % Considered channel
    channel_number = list_Nelements(kk);
    channel = rf_data(:,channel_number);
    %channel2 = channel2 ./ max(abs(channel2));
    
    %-- Measurement vector
    y = Phi*channel(:);
    param_solver.gamma = 1e-1 *norm(abs(At(y)), inf);                    % regularization parameter for the problem
    T = @(x) x;
    Tt = @(x) x;
    
    % setting different parameter for the simulation
    param_solver.verbose = 0; % display parameter
    param_solver.max_iter = 2000;       % maximum iteration
    param_solver.rel_obj = 1e-10;        % tolerance to stop iterating
    epsilon = 3e-1*norm(y);
    
    % solving the problem
    alpha_est = admm_bpcon(y, epsilon, A, At, T, Tt, param_solver);
    rf_data_rec1(:,channel_number) = Psi*alpha_est;
end

%% Beamforming in frequency
dBRange = 40;
interp = 4;

%-- Reference image
probe.pitch = probe.Pitch;
probe.c = probe.c0;
[xim, zim, bmode_compressed_ref, ~] = fkmig(rf_data, probe);
x{2} = linspace(xim(1),xim(end),length(xim)*interp);
z{2} = linspace(zim(1),zim(end),length(zim)*interp);
[X,Z] = meshgrid(x{2},z{2});
bmode_compressed_ref = interp2(xim,zim,double(bmode_compressed_ref),X,Z);

%-- Reconstructed image
[~, ~, bmode_compressed_rec1, ~] = fkmig(rf_data_rec1, probe);
bmode_compressed_rec1 = interp2(xim,zim,double(bmode_compressed_rec1),X,Z);


%-- display
figure('Color', [1 1 1])
imagesc(xim, zim, bmode_compressed_ref); colormap gray; caxis([-dBRange, 0]);
axis image;
h = colorbar;
set(h, 'YTick', [-40 -30 -20 -10 0]);
set(h, 'YTickLabel', {'-40 dB', '-30 dB', '-20 dB', '-10 dB', '0 dB'});
xlabel 'Lateral dimension [mm]'
xlabel 'Depth [mm]'
set(gca, 'FontSize',14);
export_fig('carotid_ref.pdf');
close all;

figure('Color', [1 1 1])
imagesc(xim, zim, bmode_compressed_rec1); colormap gray; caxis([-dBRange, 0]);
axis image;
xlabel 'Lateral dimension [mm]'
xlabel 'Depth [mm]'
h = colorbar;
set(h, 'YTick', [-40 -30 -20 -10 0]);
set(h, 'YTickLabel', {'-40 dB', '-30 dB', '-20 dB', '-10 dB', '0 dB'});
set(gca, 'FontSize',14);
export_fig('carotid_rec_1channel.pdf');
close all;

%% -- First channel (chosen as a reference)
channel_number = 1;
channel = rf_data(:,channel_number);

%-- Build the dictionary
pulse_pad = zeros(size(channel));
pulse_pad(1:numel(pulse),:) = pulse;
Psi = circulant(pulse_pad);

%-- Measurement matrix and measurement vector - single channel
comp_ratio = 1;
N = numel(channel);
M = round(comp_ratio*N);
Phi1 = randraw('normal', [0, 1/sqrt(M)], [M N]);
y = Phi1*channel(:);

% Setting the function
G1 = Phi1*Psi;
A =@(x) G1*x(:);
At =@(x) G1'*x(:);
param_solver.nu = eigs(G1'*G1, 1);
param_solver.gamma = 1e-2 *norm(At(y), inf);                    % regularization parameter for the problem
T = @(x) x;
Tt = @(x) x;

% setting different parameter for the simulation
param_solver.verbose = 0; % display parameter
param_solver.max_iter = 2000;       % maximum iteration
param_solver.tol = 1e-20;        % tolerance to stop iterating
epsilon = 0.01*norm(y);

% solving the problem
alpha_est1 = admm_bpcon(y, epsilon, A, At, T, Tt, param_solver);
rf_data_rec(:,channel_number) = Psi*alpha_est1;
alpha_est = alpha_est1;

%% Multichannel
% Measurement matrix for the multichannel
n_points = 100;

list_Nelements = 2:128;
for kk = 1:numel(list_Nelements)
    disp(['Reconstruction of channel ', num2str(list_Nelements(kk))]);
    
    % Considered channel
    channel_number2 = list_Nelements(kk);
    channel2 = rf_data(:,channel_number2);
    
    % Creation of the mask
    tmp  =alpha_est;
    mask = zeros(size(alpha_est));
    for ll = 1:n_points
        delta_ind = 2*round(abs(probe.xm(channel_number2) - probe.xm(channel_number2-1))*probe.fs / probe.c0);
        [~, ind] = max(abs(tmp));
        mask_tmp = zeros(numel(tmp),1);
        mask_tmp(max(ind-delta_ind, 1):min(ind+delta_ind, numel(alpha_est))) = 1;
        tmp(max(ind-delta_ind, 1):min(ind+delta_ind, numel(alpha_est))) = 0;
        mask = mask + mask_tmp;
    end
    mask(mask > 1) = 1;
    
    %-- Measurement vector
    y = Phi*channel2(:);
    
    % Setting the function
    A =@(x) G*(mask.*x(:));
    At =@(x) mask.*G'*x(:);
    param_solver.nu = pow_method(A,At,size(alpha_est),1e-4,100,0);
    param_solver.gamma = 1e-1 *norm(abs(At(y)), inf);                    % regularization parameter for the problem
    T = @(x) mask.*x;
    Tt = @(x) mask.*x;
    
    % setting different parameter for the simulation
    param_solver.verbose = 0; % display parameter
    param_solver.max_iter = 2000;       % maximum iteration
    param_solver.rel_obj = 1e-10;        % tolerance to stop iterating
    epsilon = 3e-1*norm(y);
    
    % solving the problem
    alpha_est = admm_bpcon(y, epsilon, A, At, T, Tt, param_solver);
    rf_data_rec(:,channel_number2) = Psi*alpha_est;
end

%% Beamforming in frequency
dBRange = 40;
interp = 4;

%-- Reference image
probe.pitch = probe.Pitch;
probe.c = probe.c0;
[xim, zim, bmode_compressed_ref, ~] = fkmig(rf_data, probe);
x{2} = linspace(xim(1),xim(end),length(xim)*interp);
z{2} = linspace(zim(1),zim(end),length(zim)*interp);
[X,Z] = meshgrid(x{2},z{2});
bmode_compressed_ref = interp2(xim,zim,double(bmode_compressed_ref),X,Z);

%-- Reconstructed image
[~, ~, bmode_compressed_rec, ~] = fkmig(rf_data_rec, probe);
bmode_compressed_rec = interp2(xim,zim,double(bmode_compressed_rec),X,Z);

figure('Color', [1 1 1])
imagesc(xim*1000, zim*1000, bmode_compressed_rec); colormap gray; caxis([-dBRange, 0]);
axis image;
xlabel 'Lateral dimension [mm]'
xlabel 'Depth [mm]'
h = colorbar;
set(h, 'YTick', [-40 -30 -20 -10 0]);
set(h, 'YTickLabel', {'-40 dB', '-30 dB', '-20 dB', '-10 dB', '0 dB'});
set(gca, 'FontSize',14);
export_fig('carotid_rec_multichannels.pdf');

%% Calculation of the PSNR on the normalized envelope
%-- Reference
env_ref = abs(hilbert(rf_data));
bmode_ref = env_ref / max(env_ref(:));

%-- 1 channel reconstruction
env_rec1 = abs(hilbert(rf_data_rec1));
bmode_rec1 = env_rec1 / max(env_rec1(:));

%-- multi-channel reconstruction
env_rec = abs(hilbert(rf_data_rec));
bmode_rec = env_rec / max(env_rec(:));

%-- PSNR
psnr_1channel = psnr(bmode_rec1, bmode_ref);
psnr_multi_channel = psnr(bmode_rec, bmode_ref);