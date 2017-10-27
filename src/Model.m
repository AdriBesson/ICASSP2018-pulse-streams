classdef Model < handle
    methods (Abstract)
        forward(x)
        backward(x)
    end
    methods
        function model_norm = get_norm(model, vec_size)
            tolerance = 1e-8;
            max_iter = 1e6;
            x=randn(vec_size);
            x=x/norm(x(:));           
            % Power method
            init_val = 1;
            for k = 1:max_iter
                y = model.forward(x);
                x = model.backward(y);
                model_norm = norm(x(:));
                rel_var = abs(model_norm-init_val)/init_val;
                if (rel_var < tolerance)
                    break;
                end
                init_val = model_norm;
                x = x/model_norm;
            end
        end
        function func = lsqr_func_handle(model, x, transp_flag)
            if strcmp(transp_flag, 'notransp')
                func = model.forward(x);
            else
                func = model.backward(x);
            end
        end
    end
end