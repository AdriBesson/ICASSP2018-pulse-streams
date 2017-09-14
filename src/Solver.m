classdef Solver 
    properties
        verbose
        max_iterations
        tolerance
        radius
        regularization_parameter
        init_sol
        measurement_model
        sparsity_model
        measurements
    end

    methods
        function solver = Solver(measurement_model, sparsity_model, measurements)
            solver.verbose = 0;
            solver.max_iterations = 1e3;
            solver.tolerance = 1e-10;
            solver.radius = 0;
            solver.regularization_parameter = 8e-1 *norm(measurement_model.backward(measurements), inf);
            solver.measurement_model = measurement_model;
            solver.sparsity_model = sparsity_model;
            solver.measurements = measurements;
        end
        
        function solver = set_verbose(solver, verbose)
            solver.verbose = verbose;
        end
        
        function solver = set_max_iterations(solver, max_iterations)
            solver.max_iterations = max_iterations;
        end
        
        function solver = set_tolerance(solver, tolerance)
            solver.tolerance = tolerance;
        end
        
        function solver = set_radius(solver, radius)
            solver.radius = radius;
        end
        
        function solver = set_regularization_parameter(solver, regularization_parameter)
            solver.regularization_parameter = regularization_parameter;
        end
        
        function output_coefficients = solve(solver)
            param_admm = get_admm_param(solver); 
            A = @(x) solver.measurement_model.forward(x);
            At = @(x) solver.measurement_model.backward(x);
            T = @(x) solver.sparsity_model.forward(x);
            Tt = @(x) solver.sparsity_model.backward(x);
            output_coefficients = admm_bpcon(solver.measurements, param_admm.epsilon, A, At, T, Tt, param_admm);
        end
            
    end
    methods (Access = private)
        function param_admm = get_admm_param(solver)
            param_admm.verbose = solver.verbose; 
            param_admm.max_iter = solver.max_iterations;
            param_admm.tol = solver.tolerance;
            param_admm.nu = solver.measurement_model.get_norm(size(solver.measurement_model.backward(solver.measurements)));
            param_admm.epsilon = solver.radius;
            param_admm.gamma = solver.regularization_parameter;
        end
    end
end