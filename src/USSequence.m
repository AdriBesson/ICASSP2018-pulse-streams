classdef USSequence
    properties (SetAccess = public)
        central_frequency = 5.208e6
        sampling_frequency = 20832000
        bandwidth = 0.67
        speed_of_sound = 1540
        excitation_cycle = 2.5
        impulse_response_cycle = 2
        number_elements = 128
        element_spacing = 1540/(5.208e6)
        number_time_samples = 2000
        initial_time = 0
    end

    
    methods
        
        function seq = set_central_frequency(seq, central_frequency)
            seq.central_frequency = central_frequency;
        end
        
        function seq = set_sampling_frequency(seq, sampling_frequency)
            seq.sampling_frequency = sampling_frequency;
        end
        
        function seq = set_bandwidth(seq, bandwidth)
            seq.bandwidth = bandwidth;
        end
        
        function seq = set_excitation_cycle(seq, excitation_cycle)
            seq.excitation_cycle = excitation_cycle;
        end
        
        function seq = set_impulse_response_cycle(seq, impulse_response_cycle)
            seq.impulse_response_cycle = impulse_response_cycle;
        end
        
        function seq = set_number_elements(seq, number_elements)
            seq.number_elements = number_elements;
        end
        
        function seq = set_element_spacing(seq, element_spacing)
            seq.element_spacing = element_spacing;
        end
        
        function seq = set_number_time_samples(seq, number_time_samples)
            seq.number_time_samples = number_time_samples;
        end
        
        function seq = set_initial_time(seq, initial_time)
            seq.initial_time = initial_time;
        end
        
        function pulse = estimate_received_pulse(seq)
            t0 = (-seq.impulse_response_cycle/2/seq.bandwidth/seq.central_frequency): 1/seq.sampling_frequency : (seq.impulse_response_cycle/2/seq.bandwidth/seq.central_frequency);
            te = (-seq.excitation_cycle/2/seq.central_frequency): 1/seq.sampling_frequency : (seq.excitation_cycle/2/seq.central_frequency);
            impulse_response = gauspuls(t0, seq.central_frequency, seq.bandwidth);
            excitation = square(2*pi*seq.central_frequency*te+pi/2);
            pulse = conv(conv(impulse_response,excitation),impulse_response);
            pulse = pulse / max(pulse(:));
        end
        
        function pulse_dictionary = get_pulse_dictionary(seq)
            pulse = estimate_received_pulse(seq);
            pulse_pad = zeros(seq.number_time_samples,1);
            pulse_pad(1:numel(pulse),:) = pulse;
            pulse_dictionary = circulant(pulse_pad);
        end
        
        function channel = generate_random_channel(seq, n_points)
            channel = zeros(seq.number_time_samples, 1);
            pulse = estimate_received_pulse(seq);
            
            pos_points = randi(seq.number_time_samples, [n_points, 1]);
            amp_points =  rand(size(pos_points));
            channel(pos_points) = amp_points;
            channel = conv(channel, pulse, 'same');
        end
        
        function [raw_data, points_locations_raw] = generate_rawdata(seq, points_locations, points_amplitudes, snr_awgn)
            raw_data = zeros(seq.number_time_samples, seq.number_elements);
            pulse = seq.estimate_received_pulse();
            element_locations = get_element_locations(seq);
            t_points = points_locations(:,2);
            x_points = points_locations(:,1);
            for pp = 1:numel(points_amplitudes)
                for ll = 1:seq.number_elements
                    t_Tx = t_points(pp) / seq.sampling_frequency / 2;
                    t_Rx = sqrt(1/seq.speed_of_sound^2*(element_locations(x_points(pp))-element_locations(ll))^2 + (t_points(pp) / 2 / seq.sampling_frequency)^2);
                    pos_raw_points = max(round((t_Tx + t_Rx)*seq.sampling_frequency-1), 1);
                    pos_raw_points = min(pos_raw_points, seq.number_time_samples);
                    points_locations_raw(ll,pp) = pos_raw_points;
                    raw_data( pos_raw_points ,ll) = points_amplitudes(pp);
                end
            end
            for ll = 1:seq.number_elements
                raw_data_conv(:,ll) = conv(raw_data(:,ll), pulse);
            end
            raw_data = raw_data_conv(1:seq.number_time_samples,:);
            raw_data = awgn(raw_data, snr_awgn);
        end
        
        function time_samples = get_time_samples(seq)
            time_samples = (seq.initial_time:seq.number_time_samples-1)/seq.sampling_frequency;
        end

        function element_locations = get_element_locations(seq)
            element_locations = (-seq.number_elements/2:1:seq.number_elements/2)*seq.element_spacing;
        end
        
    end
    
    
end