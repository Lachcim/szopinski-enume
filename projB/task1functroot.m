% ENUME MICHAŁ SZOPIŃSKI
% PROJECT B NUMBER 60
% TASK 1
% https://github.com/Lachcim/szopinski-enume

% define available algorithms
algorithms = {
    'bisection', @bisect;
    'Newton''s algorithm', @newton
};

% find all root brackets
interval = [2, 11];
brackets = rootbrac(@taskfunc, interval(1), interval(2));

% perform task for all available algorithms
for alg = 1:size(algorithms, 1)
    [algname, algfunc] = algorithms{alg, :};
    
    % plot function inside interval
    figure;
    grid on;
    hold on;
    title(['Approximate zeros of function (', algname, ')']);
    set(gca, 'XAxisLocation', 'origin');
    x = interval(1):0.05:interval(2);
    y = arrayfun(@taskfunc, x);
    plot(x, y);
    
    % iterate over root brackets
    for bracket = brackets    
        % find all zeros within the bracket using the given algorithm
        [zero, steps] = algfunc(@taskfunc, bracket(1), bracket(2), 1e-15);
        
        % plot steps on graph
        xline(bracket(1), 'color', [0.5, 0.5, 0.5]);
        xline(bracket(2), 'color', [0.5, 0.5, 0.5]);
        scatter(steps(1, 2:end), steps(2, 2:end), [], [0.929, 0.694, 0.125]);
        scatter(steps(1, 1), steps(2, 1), [], [0.635, 0.078, 0.184]);
        text(steps(1, 1), steps(2, 1), 'start', ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'top');
    end
    
    % finish and print graph
    hold off;
    set(gcf, 'PaperPosition', [0 0 6 4]);
    set(gcf, 'PaperSize', [6 4]);
    print(['report/', func2str(algfunc), 'zeros'], '-dpdf');
end

% the function as given in the task
function y = taskfunc(x)
    y = 0.7 * x * cos(x) - log(x + 1);
end

% finds the root brackets of a function within the given range
function brackets = rootbrac(func, rangestart, rangeend)
    % define search resolution
    resolution = (rangeend - rangestart) / 10;
    
    % start search at the start of the range
    a = rangestart;
    b = rangestart + resolution;
    brackets = double.empty(2, 0);
    
    % keep moving the interval until the range is exceeded
    while 1
        % if the function changes sign inside the interval, a bracket has been found
        if sign(func(a)) ~= sign(func(b))
            % save bracket
            brackets(:, size(brackets, 2) + 1) = [a, b];
        end
        
        % if the bracket can't be expanded, return
        if b == rangeend; return; end
        
        % check next bracket
        a = b;
        b = min(a + resolution, rangeend);
    end
end

% uses the bisection algorithm to find the root of a function within the given bracket
function [zero, steps] = bisect(func, a, b, tolerance)
    % initialize empty array of steps
    steps = double.empty(2, 0);
    
    % iterate algorithm until the error is within tolerance
    while 1
        % calculate midpoint
        zero = (a + b) / 2;
        steps(:, size(steps, 2) + 1) = [zero, func(zero)];
        
        % stop test
        if abs(func(zero)) <= tolerance; break; end
        
        % choose next sub-interval based on sign mismatch
        if sign(func(a)) ~= sign(func(zero))
            b = zero;
        else
            a = zero;
        end
    end
end

% uses Newton's algorithm to find the root of a function
function [zero, steps] = newton(func, a, b, tolerance)
    % initialize step array and calculate derivative step
    steps = double.empty(2, 0);
    step = sqrt(eps);
    
    % calculate first approximation of zero - midpoint of the bracket
    zero = (a + b) / 2;
    steps(:, size(steps, 2) + 1) = [zero, func(zero)];
    
    % iterate algorithm until the error is within tolerance
    while 1
        % calculate next approximation of zero
        derivative = (func(zero + step) - func(zero - step)) / (2 * step);
        zero = zero - func(zero) / derivative;
        steps(:, size(steps, 2) + 1) = [zero, func(zero)];
        
        % prevent divergence during approximation
        if zero < a || zero > b
            error('Divergent iteration');
        end
        
        % stop test
        if abs(func(zero)) <= tolerance; break; end
    end
end
