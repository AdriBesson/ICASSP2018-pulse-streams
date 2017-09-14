%-- Script used to display the results for SPL
clear all;
close all;
clc

%% Noiseless case - intersensor spacing 1 wavelength

%-- filnename
filename_prefix = 'resultsSPL/results_';
filnename_suffix = 'channel_synth_pulse.mat';
list_scenarios=[1, 2, 5, 10];
it = 1;

for ll=1:numel(list_scenarios)
    scenario = list_scenarios(ll);
    filename = strcat([filename_prefix, num2str(scenario), filnename_suffix]);
    load(filename);
    
    for kk = 1:size(nmse, 1)
        nmse(kk, isnan(nmse(kk,:))) = max(nmse(kk,:));
    end
    nmse_mean(ll, :) = mean(sqrt(nmse),2);
end


figure('Color', [1 1 1]);
plot(meas_ratio(1:15), nmse_mean(1, 1:15), '--o', 'LineWidth', 1.5);
hold on
plot(meas_ratio(1:15), nmse_mean(2, 1:15), '--x', 'LineWidth', 1.5);
hold on
plot(meas_ratio(1:15), nmse_mean(3, 1:15), '--d', 'LineWidth', 1.5);
hold on
plot(meas_ratio(1:15), nmse_mean(4, 1:15), '--+', 'LineWidth', 1.5);
grid on
xlabel 'Compression ratio (M/N)'
ylabel 'Normalized MSE'
legend('1channel', '2 channels', '5 channels', '10 channels');
set(gca,'fontsize',14, 'GridLineStyle', ':')
set(gcf, 'PaperPositionMode', 'auto');

%% Noiseless case - intersensor spacing 2 wavelengths

%-- Prefix and suffix of the filename
filename_prefix = 'resultsSPL/results_';
filnename_suffix = 'channel_synth_pulse_2lambda.mat';

%-- List of scenarios (number of channels involved)
list_scenarios=[1, 2, 5, 10];

%-- Loop over the scenarios
for ll=1:numel(list_scenarios)
    scenario = list_scenarios(ll);
    %-- Load the corresponding filename
    filename = strcat([filename_prefix, num2str(scenario), filnename_suffix]);
    load(filename);
    %-- Calculate the average mse over all the draws
    for kk = 1:size(nmse, 1)
        nmse(kk, isnan(nmse(kk,:))) = max(nmse(kk,:));
    end
    nmse_mean(ll, :) = mean(sqrt(nmse),2);
    
end


figure('Color', [1 1 1]);
plot(meas_ratio(1:15), nmse_mean(1,1:15), '--o', 'LineWidth', 1.5);
hold on
plot(meas_ratio(1:15), nmse_mean(2,1:15), '--x', 'LineWidth', 1.5);
hold on
plot(meas_ratio(1:15), nmse_mean(3,1:15), '--d', 'LineWidth', 1.5);
hold on
plot(meas_ratio(1:15), nmse_mean(4,1:15), '--+', 'LineWidth', 1.5);
grid on
xlabel 'Compression ratio (M/N)'
ylabel 'Normalized MSE'
legend('1channel', '2 channels', '5 channels', '10 channels');
set(gca,'fontsize',14, 'GridLineStyle', ':')
set(gcf, 'PaperPositionMode', 'auto');

%% Noisy case
%-- Single channel scenario
filename = 'resultsSPL/results_1channel_noisy.mat';
load(filename);
min_c = min(channel);
max_c = max(channel);

% Plot the results
figure('Color', [1 1 1])
plot(t*10^3, channel_noisy);
xlabel 'Time [ms]'

figure('Color', [1 1 1]);
plot(t*10^3, channel);
xlabel 'Time [ms]'
axis([0 0.07 min_c max_c])

figure('Color', [1 1 1]);
plot(t*10^3, channel_est);
axis([0 0.07 min_c max_c])
xlabel 'Time [ms]'

%-- 5 channels scenario
filename = 'resultsSPL/results_multichannel_noisy.mat';
load(filename);

% Plot the results
figure('Color', [1 1 1]);
plot(t*10^3, channel_est);
axis([0 0.07 min_c max_c])
xlabel 'Time [ms]'


%% Carotid
clear all;
filename = 'resultsSPL/results_carotid.mat';
load(filename);
dBRange = 40;

%-- Reference image
figure('Color', [1 1 1])
imagesc(xim*1000, zim*1000, bmode_compressed_ref); colormap gray; caxis([-dBRange, 0]);
axis image;
xlabel('Lateral dimension [mm]')
ylabel('Depth [mm]')
title 'Reference B-mode image'

figure('Color', [1 1 1])
imagesc(xim*1000, zim*1000, bmode_compressed_1channel_rec); colormap gray; caxis([-dBRange, 0]);
axis image;
xlabel('Lateral dimension [mm]')
ylabel('Depth [mm]')
title(['B-mode image (6% measurements)- 1 channel - PSNR= ', num2str(psnr_1channel), 'dB']);

figure('Color', [1 1 1])
imagesc(xim*1000, zim*1000, bmode_compressed_multichannel_rec); colormap gray; caxis([-dBRange, 0]);
axis image;
xlabel('Lateral dimension [mm]')
ylabel('Depth [mm]')
title(['B-mode image (6% measurements)- multichannel - PSNR= ', num2str(psnr_multichannel), 'dB']);