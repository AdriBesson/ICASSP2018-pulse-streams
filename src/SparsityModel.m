classdef SparsityModel < Model
    
    methods 

        function y = forward(sparsity_model, x)
                    y = x;
        end
        
        function x = backward(sparsity_model, y)
                    x = y;
        end
        
    end
end