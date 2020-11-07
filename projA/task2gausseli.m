% define maximum feasible equation count (10 * 2^n)
maxeqcount = 5;

% perform calculation for both sub-tasks
for task = 'ab'
    % collect error values for different equation counts
    errors = zeros(maxeqcount + 1, 2);
    
    % perform calculation for increasing number of equations
    for equationspow = 0:maxeqcount
        equationcount = 10 * 2^equationspow;
        
        % generate coefficient matrix A and vector b
        A = genmatrix(task, equationcount);
        b = genvector(task, equationcount);
        eqsys = [A, b];
        
        % perform gaussian elimination and back-substitution
        eqsys = gausseli(eqsys);
        eqsys = backsubst(eqsys);
        result = eqsys(:, size(eqsys, 2));
        
        % note the error
        error = norm(A * result - b);
        errors(equationspow + 1, :) = [equationcount, error];
    end
    
    % plot error data
    figure;
    plot(errors(:, 1), errors(:, 2), '-o');
    title(['Linear equation system solution error (subtask ', task, ')']);
    xlabel('Equation count');
    ylabel('Error');
    grid on;
end

% generates coefficient matrix for task A or B
function output = genmatrix(task, size)
    % initialize empty matrix
    output = zeros(size);
    
    if task == 'a'
        % 4's for the diagonal, 1's around it
        for i = 1:size
            output(i, i) = 4;
            if i ~= 1; output(i - 1, i) = 1; end
            if i ~= size; output(i + 1, i) = 1; end
        end
    else
        % formula for a_ij given explicitly
        for i = 1:size
            for j = 1:size
                output(i, j) = 6 / (7 * (i + j + 1));
            end
        end
    end
end

% generates the b vector for task A or B
function output = genvector(task, size)
   % initialize empty vector
   output = zeros(size, 1);
   
   if task == 'a'
       for i = 1:size
           output(i) = 4 + 0.3 * i;
       end
   else
       % 1/(3i) for even positions only
       for i = 2:2:size
           output(i) = 1 / (3 * i);
       end
   end
end

% performs gaussian elimnation
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
