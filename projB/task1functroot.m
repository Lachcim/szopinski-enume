% ENUME MICHAŁ SZOPIŃSKI
% PROJECT B NUMBER 60
% TASK 1
% https://github.com/Lachcim/szopinski-enume

for bracket = rootbrac(@taskfunc, 2, 11)
    zero = bisect(@taskfunc, bracket(1), bracket(2));
    disp(zero);
    zero2 = newton(@taskfunc, bracket(1), bracket(2));
    disp(zero2);
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
function zero = bisect(func, a, b)
    % iterate algorithm until the result range decreases below the threshold
    while (b - a) > 1e-12
        % calculate midpoint
        midpoint = (a + b) / 2;
        
        % choose next sub-interval based on sign mismatch
        if sign(func(a)) ~= sign(func(midpoint))
            b = midpoint;
        else
            a = midpoint;
        end
    end
    
    % pick midpoint of the final range as the result
    zero = (a + b) / 2;
end

% uses Newton's algorithm to find the root of a function
function zero = newton(func, a, b)
    % calculate square root of epsilon for derivative calculation
    step = sqrt(eps);
    
    % calculate first approximation of zero - midpoint of the bracket
    zero = (a + b) / 2;
    
    for n = 1:55
        % calculate next approximation of zero
        derivative = (func(zero + step) - func(zero - step)) / (2 * step);
        zero = zero - func(zero) / derivative;
        
        % prevent divergence during approximation
        if zero < a || zero > b
            error('Divergent iteration');
        end
    end
end
