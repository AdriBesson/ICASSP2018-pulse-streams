%-- Script used to display the results for SPL
clear all;
close all;
clc
%% Noiseless case
%-- Single channel scenario
filename = 'resultsSPL/results_1channel_synth_pulse.mat';
load(filename);

for kk = 1:size(nmse, 1)
    nmse(kk, isnan(nmse(kk,:))) = max(nmse(kk,:));
end
nmse_mean_1 = mean(sqrt(nmse),2);

%-- 2 channels scenario
filename = 'resultsSPL/results_2channels_synth_pulse.mat';
load(filename);

for kk = 1:size(nmse, 1)
    nmse(kk, isnan(nmse(kk,:))) = max(nmse(kk,:));
end
nmse_mean_2 = mean(sqrt(nmse),2);

%-- 5 channels scenario
filename = 'resultsSPL/results_5channels_synth_pulse.mat';
load(filename);

for kk = 1:size(nmse, 1)
    nmse(kk, isnan(nmse(kk,:))) = max(nmse(kk,:));
end
nmse_mean_5 = mean(sqrt(nmse),2);

%-- 10 channels scenario
filename = 'resultsSPL/results_10channels_synth_pulse.mat';
load(filename);

for kk = 1:size(nmse, 1)
    nmse(kk, isnan(nmse(kk,:))) = max(nmse(kk,:));
end
nmse_mean_10 = mean(sqrt(nmse),2);


figure('Color', [1 1 1]);
plot(meas_ratio, nmse_mean_1, '--o', 'LineWidth', 1.5);
hold on 
plot(meas_ratio, nmse_mean_2, '--x', 'LineWidth', 1.5);
hold on 
plot(meas_ratio, nmse_mean_5, '--d', 'LineWidth', 1.5);
hold on 
plot(meas_ratio, nmse_mean_10, '--+', 'LineWidth', 1.5);
grid on
xlabel 'Compression ratio (M/N)'
ylabel 'Normalized MSE'
legend('1channel', '2 channels', '5 channels', '10 channels');
set(gca,'fontsize',14, 'GridLineStyle', ':')
set(gcf, 'PaperPositionMode', 'auto');

%% Noiseless case 2 lambda
%-- Single channel scenario
filename = 'resultsSPL/results_1channel_synth_pulse.mat';
load(filename);

for kk = 1:size(nmse, 1)
    nmse(kk, isnan(nmse(kk,:))) = max(nmse(kk,:));
end
nmse_mean_1 = mean(sqrt(nmse),2);

%-- 2 channels scenario
filename = 'resultsSPL/results_2channels_synth_pulse_2lambda.mat';
load(filename);

for kk = 1:size(nmse, 1)
    nmse(kk, isnan(nmse(kk,:))) = max(nmse(kk,:));
end
nmse_mean_2 = mean(sqrt(nmse),2);

%-- 5 channels scenario
filename = 'resultsSPL/results_5channels_synth_pulse_2lambda.mat';
load(filename);

for kk = 1:size(nmse, 1)
    nmse(kk, isnan(nmse(kk,:))) = max(nmse(kk,:));
end
nmse_mean_5 = mean(sqrt(nmse),2);

%-- 10 channels scenario
filename = 'resultsSPL/results_10channels_synth_pulse.mat';
load(filename);

for kk = 1:size(nmse, 1)
    nmse(kk, isnan(nmse(kk,:))) = max(nmse(kk,:));
end
nmse_mean_10 = mean(sqrt(nmse),2);


figure('Color', [1 1 1]);
plot(meas_ratio, nmse_mean_1, '--o', 'LineWidth', 1.5);
hold on 
plot(meas_ratio, nmse_mean_2, '--x', 'LineWidth', 1.5);
hold on 
plot(meas_ratio, nmse_mean_5, '--d', 'LineWidth', 1.5);
hold on 
plot(meas_ratio, nmse_mean_10, '--+', 'LineWidth', 1.5);
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


