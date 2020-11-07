% perform calculation for both sub-tasks
for task = 'ab'
    % generate coefficient matrix A and vector b
    A = genmatrix(task, 5);
    b = genvector(task, 5);
    eqsys = [A, b];
    
    disp(eqsys)
    eqsys = gausselim(eqsys);
    disp(eqsys)
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
function eqsys = gausselim(eqsys)
    % for every column, eliminate (size - column) coefficients
    for col = 1:size(eqsys, 1)
        for row = (col + 1):size(eqsys, 1)
            % calculate reductor constant and perform reduction
            reductor = eqsys(row, col) / eqsys(col, col);
            eqsys(row, :) = eqsys(row, :) - eqsys(col, :) * reductor;
        end
    end
end
