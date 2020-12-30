% ENUME MICHAŁ SZOPIŃSKI
% PROJECT C NUMBER 60
% TASK 1
% https://github.com/Lachcim/szopinski-enume

% use functions from project A if not present in the working directory
if ~exist('qrdecomp', 'var')
    addpath('../projA');
end

% define function data points
taskfunc = (-5:5)';
taskfunc(:, 2) = [
    -15.2991;
    -11.9874;
    -7.8757;
    -5.7178;
    -3.3653;
    -2.5691;
    -3.3150;
    -6.2274;
    -10.7044;
    -19.1618;
    -30.7795
];

% find the approximating polynomial of the given degree
function factors = approximate(func, polydeg)
    % define the A matrix used for calculating Gram's matrix
    A = zeros(size(func, 1), polydeg + 1);
    
    % calculate cells of A using natural basis
    for row = 1:size(A, 1)
        for col = 1:size(A, 2)            
            A(row, col) = func(row, 1) ^ (col - 1);
        end
    end
    
    % solve least-square problem using QR decomposition
    [Q, R] = qrdecomp(A);
    eqsys = R(1:size(R, 2), :);
    eqsys(:, end + 1) = Q' * func(:, 2);
    factors = backsubst(eqsys);
end
