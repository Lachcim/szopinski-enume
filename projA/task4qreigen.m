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
disp(size(errors, 2));
[eigenvalues, errors] = eigenshifts(A);
disp(eigenvalues);
disp(size(errors, 2));

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
        errors(size(errors, 2) + 1) = maxnonzero;
        if (maxnonzero <= 1e-6); break; end
    end
    
    % convert eigenvalue matrix to vector
    eigenvalues = diag(A)';
end

% finds the eigenvalues of a matrix using QR with shifts
function [eigenvalues, errors] = eigenshifts(A)
    % initialize empty output
    eigenvalues = double.empty(1, 0);
    errors = double.empty(1, 0);
    
    % consider increasingly smaller sub-matrices
    matsize = size(A, 1);
    while matsize >= 2
        % find one eigenvalue
        while 1
            % find eigenvalue of the lower right corner
            corner = A((matsize - 1):matsize, (matsize - 1):matsize);
            cornereigen = eigenoftwo(corner);
            
            % shift and iterate algorithm
            A = A - eye(matsize) * cornereigen;
            [Q, R] = qrdecomp(A);
            A = R * Q + eye(matsize) * cornereigen;
            
            % move on once the row has been zeroed out
            maxnonzero = max(A(matsize, 1:(matsize - 1)));
            errors(size(errors, 2) + 1) = maxnonzero;
            
            if (maxnonzero <= 1e-6)
                % remember discovered eigenvalue
                eigenvalues(size(eigenvalues, 2) + 1) = A(matsize, matsize);
                break;
            end
        end
        
        % deflate matrix
        matsize = matsize - 1;
        A = A(1:matsize, 1:matsize);
    end
    
    % include final eigenvalue
    eigenvalues(size(eigenvalues, 2) + 1) = A(1, 1);
end

% finds the eigenvalue of a 2x2 matrix that is closer to the lower right corner
function eigen = eigenoftwo(A)
    % solve characteristic equation of matrix to find its eigenvalues
    delta = (A(1) + A(4)) ^ 2 - 4 * (A(1) * A(4) - A(2) * A(3));
    sqrtdelta = sqrt(delta);
    eigen1 = ((A(1) + A(4)) - sqrtdelta) / 2;
    eigen2 = ((A(1) + A(4)) + sqrtdelta) / 2;
    
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
