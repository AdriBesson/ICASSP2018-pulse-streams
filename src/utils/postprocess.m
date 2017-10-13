function [ bmode ] = postprocess( rf_data )
% Envelope detection
env  = double(abs(hilbert(double(rf_data))));

% Normalization
bmode  = 100*env./max(abs(env(:)));

end

