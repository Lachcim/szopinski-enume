% ENUME MICHAŁ SZOPIŃSKI
% PROJECT A NUMBER 62
% TASK 2
% https://github.com/Lachcim/szopinski-enume

% define maximum feasible equation count (10 * 2^n)
maxeqcount = 5;

% perform calculation for both sub-tasks
for task = 'ab'
    % collect error values for different equation counts
    errors = zeros(maxeqcount + 1, 2);
    
    % perform calculation for increasing number of equations
    for equationspow = 0:maxeqcount
        equationcount = 10 * 2 ^ equationspow;
        
        % generate coefficient matrix A and vector b
        A = genmatrix(task, equationcount);
        b = genvector(task, equationcount);
        eqsys = [A, b];
        
        % perform gaussian elimination and back-substitution
        eqsys = gausseli(eqsys);
        eqsys = backsubst(eqsys);
        result = eqsys(:, size(eqsys, 2));
        
        % note the error
        errorvector = A * result - b;
        error = norm(errorvector);
        errors(equationspow + 1, :) = [equationcount, error];
        
        % additional calculations for 10 equations
        if equationspow == 0
            % print solution
            disp(['Solution to subtask ', task ', 10 equations:']);
            disp([(1:10)', result, errorvector]);
            disp(['Error: ', num2str(error)]);
            
            % residual correction: calculate deltax and correct result
            for i = 1:10
                eqsys = [A, A * result - b];
                eqsys = gausseli(eqsys);
                eqsys = backsubst(eqsys);
                
                deltax = eqsys(:, size(eqsys, 2));
                result = result - deltax;
            end
            
            % calculate error of corrected result
            errorvector = A * result - b;
            error = norm(errorvector);

            % print corrected result
            disp('With residual correction:');
            disp([(1:10)', result, errorvector]);
            disp(['Error: ', num2str(error)]);
        end
    end
    
    % plot error data
    figure;
    plot(errors(:, 1), errors(:, 2), '-o');
    title(['Linear equation system solution error (subtask ', task, ')']);
    xlabel('Equation count');
    ylabel('Error');
    grid on;
    set(gcf, 'PaperPosition', [0 0 6 4]);
    set(gcf, 'PaperSize', [6 4]);
    print(['report/', task, 'error'], '-dpdf')
end

% performs Gaussian elimination
function eqsys = gausseli(eqsys)
    % for every column, eliminate (size - column) coefficients
    for col = 1:size(eqsys, 1)
        % partial pivoting - find best row
        bestrow = 0;
        bestscore = 0;
        for row = col:size(eqsys, 1)
            score = abs(eqsys(row, col));
            if score > bestscore
                bestrow = row;
                bestscore = score;
            end
        end
        
        % swap first row and best row
        eqsys([col, bestrow], :) = eqsys([bestrow, col], :);
        
        for row = (col + 1):size(eqsys, 1)
            % calculate reductor constant and perform reduction
            reductor = eqsys(row, col) / eqsys(col, col);
            eqsys(row, :) = eqsys(row, :) - eqsys(col, :) * reductor;
            
            % simulate perfect reduction
            eqsys(row, col) = 0;
        end
    end
end

% reduces triangular matrix to diagonal matrix
function eqsys = backsubst(eqsys)
    for col = size(eqsys, 1):-1:1
        % normalize diagonal coefficients to 1
        eqsys(col, :) = eqsys(col, :) / eqsys(col, col);
        
        % eliminate factor from other rows
        for row = (col - 1):-1:1
            reductor = eqsys(row, col) / eqsys(col, col);
            eqsys(row, :) = eqsys(row, :) - eqsys(col, :) * reductor;
        end
    end
end
