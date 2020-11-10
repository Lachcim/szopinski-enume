% ENUME MICHAŁ SZOPIŃSKI
% PROJECT A NUMBER 62
% TASK 4
% https://github.com/Lachcim/szopinski-enume

A = [ 1,  1, 9;
      2, -1, 0;
     -2,  4, 0];
[Q, R] = qrdecomp(A);
disp(Q);
disp(R);

function [Q, R] = qrdecomp(A)
    % initialize empty matrices
    Q = zeros(size(A));
    R = eye(size(A, 2));
    
    % calculate Q and R using Gram-Schmidt algorithm
    for col = 1:size(A, 2)
        % initialize column as its unorthogonalized counterpart
        Q(:, col) = A(:, col);
        
        % iterate over previous columns to obtain the orthogonalized column
        for prev = 1:(col - 1)
            % calculate the R cell for this column pair
            R(prev, col) = dot(Q(:, prev), A(:, col)) / dot(Q(:, prev), Q(:, prev));
            
            % from the current column, subtract the previous column times the R cell
            Q(:, col) = Q(:, col) - R(prev, col) * Q(:, prev);
        end
    end
end
