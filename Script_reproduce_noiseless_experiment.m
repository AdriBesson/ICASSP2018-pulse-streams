% Regenerate the noiseless experiment
addpath(genpath('src'));
cd('Noiseless recovery');

%1 channel
run('Script_1channel_synthetic_noiseless_class.m')

% Multichannel
run('Script_multichannel_synthetic_noiseless_class.m')