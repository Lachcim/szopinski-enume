% ENUME MICHAŁ SZOPIŃSKI
% PROJECT B NUMBER 60
% TASK 2
% https://github.com/Lachcim/szopinski-enume

% find all root brackets
interval = [1, 7];
brackets = rootbrac(@polynomial, interval(1), interval(2));
[root, steps] = mm1(@polynomial, brackets(1, 1), brackets(2, 1), 1e-15);
disp(steps);

% find roots of polynomial using MM1
function [zero, steps] = mm1(func, a, b, tolerance)
    % initialize output
    steps = double.empty(2, 0);
    
    % define the three approximation points
    apprx = [a, b, (a + b) / 2];
    apprxval = arrayfun(func, apprx);
    
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
