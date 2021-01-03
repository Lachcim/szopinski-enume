% ENUME MICHAŁ SZOPIŃSKI
% PROJECT C NUMBER 60
% TASK 2
% https://github.com/Lachcim/szopinski-enume

% functions of the ODE system and initial values
sysfuncts = {
    @(x) x(2) + x(1) * (0.5 - x(1)^2 - x(2)^2);
    @(x) -x(1) + x(2) * (0.5 - x(1)^2 - x(2)^2)
};
initvalues = [0; 0];

eulerode(sysfuncts, initvalues, 0, 15, 1)

function x = eulerode(functs, init, a, b, stepsize)
    % set initial values as start points of output
    x = init;
    
    % build output based on preceding values
    stepcount = ceil((b - a)/stepsize);
    for step = 1:stepcount
        for eqnum = 1:size(functs, 1)
            prevderiv = functs{eqnum}(x(:, step));
            
            x(eqnum, step + 1) = x(eqnum, step) + stepsize * prevderiv;
        end
    end
    
    % append arguments to output
    x = [a:stepsize:(stepcount * stepsize); x];
end
