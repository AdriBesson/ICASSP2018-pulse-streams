function support = get_support(us_seq, alpha, size_support)
    delta= round(abs(us_seq.element_spacing)*us_seq.sampling_frequency / us_seq.speed_of_sound);
    support = zeros(size_support, 1);
    for kk = 1:size_support
        [~, ind] = max(abs(alpha));
        alpha(max(ind-delta, 1):min(ind+delta, numel(alpha))) = 0;
        support(kk) = ind;
    end
end