% ENUME MICHAŁ SZOPIŃSKI
% PROJECT A NUMBER 62
% https://github.com/Lachcim/szopinski-enume

% performs QR decomposition of a matrix
function [Q, R] = qrdecomp(A)
    % initialize empty matrices
    Q = zeros(size(A));
    R = eye(size(A, 2));
    
    % modified Gram-Schmidt, use each column to orthogonalize the ones in front of it
    for col = 1:size(A, 2)
        % by the time we've reached this column, it's already been orthogonalized
        Q(:, col) = A(:, col);
        
        % calculate current column dot product for R
        coldot = dot(Q(:, col), Q(:, col));
        
        % orthogonalize further columns
        for next = (col + 1):size(A, 2)
            % calculate R cell for this column pair
            R(col, next) = dot(Q(:, col), A(:, next)) / coldot;
            
            % orthogonalize column
            A(:, next) = A(:, next) - R(col, next) * Q(:, col);
        end
    end
    
    % normalize matrix
    normalizer = zeros(size(Q));
    for col = 1:size(Q, 2)
        colnorm = norm(Q(:, col));
        normalizer(col, col) = colnorm;
        Q(:, col) = Q(:, col) / colnorm;
    end
    R = normalizer * R;
end
