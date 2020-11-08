% ENUME MICHAŁ SZOPIŃSKI
% PROJECT A NUMBER 62
% TASK 3
% https://github.com/Lachcim/szopinski-enume

% splits a matrix a into lower, upper and diagonal matrices
function [lower, upper, diagonal] = splitmatrix(input)
    % allocate empty matrices
    lower = zeros(size(input));
    upper = zeros(size(input));
    diagonal = zeros(size(input));
    
    % copy each element to the right output matrix
    for i = 1:size(input, 1)
        for j = 1:size(input, 2)
            if i > j; lower(i, j) = input(i, j);
            elseif i < j; upper(i, j) = input(i, j);
            else; diagonal(i, j) = input(i, j);
            end
        end
    end
end
