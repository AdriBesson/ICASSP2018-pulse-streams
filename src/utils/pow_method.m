function val = pow_method(A, At, im_size, tol, max_iter, verbose)
% Power method used to compute the maximum eigenvalue of AtA
%
% val = pow_method(A, At, im_size, tol, max_iter, verbose)
%
%
% The input argument param contains the following fields:
%
%   - A: Operator
%
%   - At: Adjoint of A
%
%	- im_size: size of the input of A 
%
%	- tol: tolerance of the algorithm (pow method stops when relative error < tol)
%
%   - max_iter: maximum number of iterations of the power method
%
%   - verbose: 0 (no diplay) / 1 (display)
%

% Initialization   
x=randn(im_size);
x=x/norm(x(:));
init_val=1;

% Loop of the power method
for k=1:max_iter
    y=A(x);
    x=At(y);
    val=norm(x(:));
    rel_var=abs(val-init_val)/init_val;
    if (verbose > 0)
        fprintf('Iter = %i, norm = %e \n',k,val);
    end
    if (rel_var < tol)
        break;
    end
    init_val=val;
    x=x/val;
    
end


end

