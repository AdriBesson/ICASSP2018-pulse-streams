% Regenerate the noiseless experiment
cd('Noiseless recovery');

%1 channel
run('Script_1channel_synthetic_noiseless.m')

% Multichannel
run('Script_multichannel_synthetic_noiseless.m')