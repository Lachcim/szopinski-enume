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
interval = [0; 15];

% define available algorithms
algorithms = {
    'RK4', @rk4, [0.01, 0.013408];
};

% solve ODE using different algorithms and step sizes
for alg = 1:size(algorithms, 1)
    [algname, algfunc, stepsizes] = algorithms{alg, :};
    
    % solve using the given algorithm for each step size
    stepresults = cell(size(stepsizes, 2), 3);
    stepnames = {'optimal step', 'larger step'};
    for stepno = 1:size(stepsizes, 2)
        result = algfunc(sysfuncts, initvalues, interval(1), interval(2), stepsizes(stepno));
        stepresults(stepno, :) = { ...
            stepsizes(stepno), ...
            stepnames{stepno}, ...
            result};
    end
    
    disp(stepresults);
end
