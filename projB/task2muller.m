% ENUME MICHAŁ SZOPIŃSKI
% PROJECT B NUMBER 60
% TASK 2
% https://github.com/Lachcim/szopinski-enume

% define available algorithms
algorithms = {
    'MM1', @mm1;
    'MM2', @mm2
};

% find all real root brackets
interval = [1, 7];
brackets = rootbrac(@polynomial, interval(1), interval(2));

% find and graph real roots using both algorithms
printroots(@polynomial, algorithms, interval, brackets, ...
    'Approximate real roots of polynomial', 'realroots');
printcomplex(@polynomial, algorithms, [-1+i, 0], ...
    'Approximate complex roots of polynomial', 'complexroots');

% find roots of polynomial using MM1
function [zero, steps] = mm1(func, a, b, tolerance)
    % define the three approximation points
    apprx = [a, b, (a + b) / 2];
    apprxval = arrayfun(func, apprx);
    
    % initialize output
    steps = [apprx(3); func(apprx(3))];
    
    % iterate algorithm until the error is within tolerance
    while abs(apprx(3) - apprx(2)) > tolerance
        % prepare linear equation system to find parabola
        z0 = apprx(1) - apprx(3);
        z1 = apprx(2) - apprx(3);
        diff0 = apprxval(1) - apprxval(3);
        diff1 = apprxval(2) - apprxval(3);

        % solve equation system
        eqsysleft = [z0 ^ 2, z0; z1 ^ 2, z1];
        eqsysright = [diff0; diff1];
        eqsyssol = eqsysleft \ eqsysright;

        % define approximation parabola
        a = eqsyssol(1);
        b = eqsyssol(2);
        c = apprxval(3);

        % find roots of parabola
        zplus = -2 * c / (b + sqrt(b ^ 2 - 4 * a * c));
        zminus = -2 * c / (b - sqrt(b ^ 2 - 4 * a * c));

        % choose root closer to current approximation
        if abs(zplus) < abs(zminus)
            newapprx = apprx(3) + zplus;
        else
            newapprx = apprx(3) + zminus;
        end
        
        % update answer
        zero = newapprx;
        steps(:, size(steps, 2) + 1) = [zero, func(zero)];
        
        % eliminate the most distant of the three approximations
        worstapprxindex = -1;
        worstapprxdiff = 0;
        for i = 1:size(apprx, 2)
            diff = abs(apprx(i) - newapprx);
            if diff > worstapprxdiff
                worstapprxindex = i;
                worstapprxdiff = diff;
            end
        end
        
        % delete old approximation and append new one
        apprx(worstapprxindex) = [];
        apprx(3) = newapprx;
        apprxval = arrayfun(func, apprx);
    end
end

% find roots of polynomial using MM2
function [approx, steps] = mm2(func, a, b, tolerance)
    % define current and (dummy) previous approximation point
    approx = (a + b) / 2;
    prevapprox = approx + b - a;
    
    % initialize output
    steps = [approx; func(approx)];
    
    % iterate algorithm until the error is within tolerance
    % the error is defined as the diff between the prev and the current approx
    while abs(approx - prevapprox) > tolerance
        % calculate approximating parabola using first and second derivative
        c = func(approx);
        b = deriv(func, approx, 1);
        a = deriv(func, approx, 2) / 2;
        
        % find roots of parabola
        zplus = -2 * c / (b + sqrt(b ^ 2 - 4 * a * c));
        zminus = -2 * c / (b - sqrt(b ^ 2 - 4 * a * c));
        
        % choose root closer to current approximation
        if abs(zplus) < abs(zminus)
            newapprox = approx + zplus;
        else
            newapprox = approx + zminus;
        end
        
        % update answer and prev approx
        prevapprox = approx;
        approx = newapprox;
        steps(:, size(steps, 2) + 1) = [approx, func(approx)];
    end
end

% calculate the nth derivative of func at x
function y = deriv(func, x, deg)
    % base case: zeroth derivative
    if deg == 0
        y = func(x);
        return
    end
    
    % recurse to find the nth derivative
    step = sqrt(eps);
    y = (deriv(func, x + step, deg - 1) - deriv(func, x - step, deg - 1)) / (2 * step);
end
