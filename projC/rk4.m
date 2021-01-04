% ENUME MICHAŁ SZOPIŃSKI
% PROJECT C NUMBER 60
% https://github.com/Lachcim/szopinski-enume

% solve ODE system using RK4 with constant step size
function x = rk4(functs, init, interval, stepsize)
    % set initial values as start points of output
    x = init;
    
    % build output based on preceding values
    stepcount = ceil((interval(2) - interval(1)) / stepsize);
    for step = 1:stepcount
        % obtain the preceding function values
        stepval = x(:, step);
        
        for eqnum = 1:size(functs, 1)
            % generic single-step iteration
            phi = rk4phi(functs{eqnum}, stepval, stepsize);
            x(eqnum, step + 1) = x(eqnum, step) + stepsize * phi;
        end
    end
    
    % append arguments to output
    x = [interval(1):stepsize:(stepcount * stepsize); x];
end
