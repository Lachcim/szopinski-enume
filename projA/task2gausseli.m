q = genvector('2', 10);

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
