% ENUME MICHAŁ SZOPIŃSKI
% PROJECT A NUMBER 62
% TASK 3
% https://github.com/Lachcim/szopinski-enume

% define matrix A and vector b as specified in the task
taskA =	[	18,		2,		-3,		1;
			2,		-25,	5,		-18;
			1,		3,		13,		-8;
			1,		1,		-2,		-10;	];
taskb = [7, 12, 24, 20]';

disp(taskA \ taskb);
disp(jacobi(taskA, taskb));

function x = jacobi(A, b)
	% split input matrix and create step zero result vector
	[lower, upper, invdiagonal] = splitmatrix(A);
	x = ones(size(A, 1), 1);
	
	M = -invdiagonal * (lower + upper);
	w = invdiagonal * b;
	
	% execute the algorithm until the desired accuracy is achieved
	while 1
		% calculate next iteration
		x = M * x + w;
		
		% calculate error
		errorvector = A * x - b;
        error = norm(errorvector);
		
        % stop iteration when the error drops below the threshold
		if error < 1e-9; break; end
	end
end

% splits a matrix a into lower, upper and inverse diagonal matrices
function [lower, upper, invdiagonal] = splitmatrix(input)
    % allocate empty matrices
    lower = zeros(size(input));
    upper = zeros(size(input));
    invdiagonal = zeros(size(input));
    
    % copy each element to the right output matrix
    for i = 1:size(input, 1)
        for j = 1:size(input, 2)
            if i > j; lower(i, j) = input(i, j);
            elseif i < j; upper(i, j) = input(i, j);
            else; invdiagonal(i, j) = 1 / input(i, j);
            end
        end
    end
end
