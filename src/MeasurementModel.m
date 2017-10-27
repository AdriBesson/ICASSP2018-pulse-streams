classdef MeasurementModel < Model
    properties 
        measurement_matrix
        mask = 1
    end
    
    methods 
        function measurement_model = MeasurementModel(measurement_matrix)
            measurement_model.measurement_matrix = measurement_matrix;
        end
        
        function measurement_model = set_mask(measurement_model, mask)
            measurement_model.mask = mask;
        end
        
        function y = forward(measurement_model, x)
            y = measurement_model.measurement_matrix*((measurement_model.mask).*x(:));
        end
        
        function x = backward(measurement_model, y)
            x = (measurement_model.mask).*(transpose(measurement_model.measurement_matrix)*y(:));
        end
        
    end
end