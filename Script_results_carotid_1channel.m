%-- Script used to reconstruct a real carotid
clear all;
close all;
clc
addpath(genpath('src'));
addpath(genpath('data'));

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
M = round(0.08*N);
Phi = randraw('normal', [0, 1/sqrt(M)], [M N]);
G = Phi*Psi;
A =@(x) G*x(:);
At =@(x) G'*x(:);
param_solver.nu = eigs(G'*G,1);
list_Nelements = 1:128;

%-- Reconstructed raw data
rf_data_rec = zeros(size(rf_data));

for kk = 1:numel(list_Nelements)
    disp(['Reconstruction of channel ', num2str(list_Nelements(kk))]);
    
    % Considered channel
    channel_number = list_Nelements(kk);
    channel = rf_data(:,channel_number);
    %channel2 = channel2 ./ max(abs(channel2));
    
    %-- Measurement vector
    y = Phi*channel(:);
    param_solver.gamma = 5e-1 *norm(abs(At(y)), inf);                    % regularization parameter for the problem
    T = @(x) x;
    Tt = @(x) x;
    
    % setting different parameter for the simulation
    param_solver.verbose = 0; % display parameter
    param_solver.max_iter = 1000;       % maximum iteration
    param_solver.rel_obj = 1e-10;        % tolerance to stop iterating
    epsilon = 3e-1*norm(y);
    
    % solving the problem
    alpha_est = admm_bpcon(y, epsilon, A, At, T, Tt, param_solver);
    rf_data_rec(:,channel_number) = Psi*alpha_est;
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


%-- display
figure('Color', [1 1 1])
imagesc(xim, zim, bmode_compressed_ref); colormap gray; caxis([-dBRange, 0]);
axis image;

figure('Color', [1 1 1])
imagesc(xim, zim, bmode_compressed_rec); colormap gray; caxis([-dBRange, 0]);
axis image;