% ENUME MICHAŁ SZOPIŃSKI
% PROJECT A NUMBER 62
% TASK 4
% https://github.com/Lachcim/szopinski-enume

A = [1, 1, 7, 5, 2;
     1, 8, 5, 4, 4;
     7, 5, 0, 8, 8;
     5, 4, 8, 0, 8;
     2, 4, 8, 8, 1];

[eigenvalues, errors] = eigennoshifts(A);
disp(eigenvalues);

% finds the eigenvalues of a matrix using the QR method without shifts
function [eigenvalues, errors] = eigennoshifts(A)
    % initialize empty error array
    errors = double.empty(1, 0);
    
    while 1
        % converge to eigenvalue diagonal matrix
        [Q, R] = qrdecomp(A);
        A = R * Q;
        
        % iterate until all non-diagonal elements are below the threshold
        nondiag = A - diag(diag(A));
        maxnonzero = max(max(nondiag));
        
        % log error
        errors(size(errors, 2) + 1) = maxnonzero;
        
        % stop iteration
        if (maxnonzero <= 1e-6); break; end
    end
    
    % convert eigenvalue matrix to vector
    eigenvalues = diag(A);
end

% finds the eigenvalue of a 2x2 matrix that is closer to the lower right corner
function eigen = eigenoftwo(A)
    % solve characteristic equation of matrix to find its eigenvalues
    delta = (A(1) + A(4)) ^ 2 - 4 * (A(1) * A(4) - A(2) * A(3));
    eigen1 = ((A(1) + A(4)) - sqrt(delta)) / 2;
    eigen2 = ((A(1) + A(4)) + sqrt(delta)) / 2;
    
    % return value closer to lower right corner
    if abs(A(4) - eigen1) < abs(A(4) - eigen2)
        eigen = eigen1;
    else
        eigen = eigen2;
    end
end

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
