function [xsol] = mb_iht(y, A, At, mask, s, tol, verbose)
% Solve model based iterative hard thresholding algorithm

% Length of the original signal
N = length(At(y));

% length of the measurements
n = length(y);

% Initial solution
xsol = zeros(N,1);

% Initial residue 
r = y;

%Main loop. Sequential.
for t = 1:param.max_iter
    % Pre threshold value
    gamma = xsol + At(r);
    
    % Find the s-th largest value of gamma
    gamma_sorted = sort(abs(mask.*gamma), 'descend');
    threshold = gamma_sorted(s);
    
    % Hard threshold of the masked signal
    xhat_prev = xhat;
    xhat(mask) = wthresh(mask.*gamma, 'h', threshold);
    
    %Check relative change of objective function   
    rel_fval = abs(xhat - xhat_prev)/xhat;
    
    % Update the residuals
    r = y - A(xhat);
    
    %Log
    if (verbose >= 1)
        fprintf('Iter %i\n',t);
        fprintf('relative residual = %e\n\n', rel_fval);
    end
    
    
    %Global stopping criteria
    if (rel_fval < tol)
        break;
    end
end
end

