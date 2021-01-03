% ENUME MICHAŁ SZOPIŃSKI
% PROJECT C NUMBER 60
% TASK 2
% https://github.com/Lachcim/szopinski-enume

% functions of the ODE system and initial values
sysfuncts = {
    @(x) x(2) + x(1) * (0.5 - x(1)^2 - x(2)^2);
    @(x) -x(1) + x(2) * (0.5 - x(1)^2 - x(2)^2)
};
initvalues = [0; 12];

rk4(sysfuncts, initvalues, 0, 15, 0.013408);

function x = rk4(functs, init, a, b, stepsize)
    % set initial values as start points of output
    x = init;
    
    % build output based on preceding values
    stepcount = ceil((b - a)/stepsize);
    for step = 1:stepcount
        % obtain the preceding function values
        stepval = x(:, step);
        
        for eqnum = 1:size(functs, 1)
            % obtain the function for this equation
            fun = functs{eqnum};
            
            % calculate the phi for RK4
            k1 = fun(stepval);
            k2 = fun(stepval + 0.5 * stepsize * k1);
            k3 = fun(stepval + 0.5 * stepsize * k2);
            k4 = fun(stepval + stepsize * k3);
            phi = (k1 + 2 * k2 + 2 * k3 + k4) / 6;
            
            % generic single-step iteration
            x(eqnum, step + 1) = x(eqnum, step) + stepsize * phi;
        end
    end
    
    % append arguments to output
    x = [a:stepsize:(stepcount * stepsize); x];
end
