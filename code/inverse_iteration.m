function c = inverse_iteration(alpha, beta, lambda_i, delta)
% INVERSE_ITERATION  Finds the eigenvector of Gamma for lambda_i
% using inverse iteration with Thomas algorithm
% Following Li et al. (2000), Section 4
%
% INPUT:
%   alpha    - diagonal of Gamma (m1 x 1)
%   beta     - subdiagonal of Gamma (m1-1 x 1)
%   lambda_i - known eigenvalue (from bisection_eig)
%   delta    - perturbation to avoid singularity (default 1e-10)
%
% OUTPUT:
%   c - normalized eigenvector (m1 x 1)

    if nargin < 4
        delta = 1e-10;
    end

    m1 = length(alpha);

    %  Shifted matrix: A = Gamma - lambda_i*I + delta*I 
    % Shifted diagonal
    d = alpha - lambda_i + delta;  % (m1 x 1)
    % Subdiagonal remains beta

    %  Random normalized initial vector 
    c = randn(m1, 1);
    c = c / norm(c);

    %  Inverse iteration: 3 iterations are always sufficient 
    % if lambda_i is known with precision 1e-10
    max_iter = 3;

    for iter = 1:max_iter

        % Solve (Gamma - lambda_i*I + delta*I) * c_new = c
        % i.e.: A * c_new = c
        % with A tridiagonal symmetric → Thomas algorithm

        c_new = thomas_tridiag(d, beta, c);

        % Normalize
        c_new_norm = norm(c_new);

        if c_new_norm < 1e-14
            % Numerically zero — convergence problem
            warning('inverse_iteration: nearly zero norm at iter %d', iter);
            break
        end

        c = c_new / c_new_norm;

    end

end


function x = thomas_tridiag(d, beta, b)
% THOMAS_TRIDIAG  Solves tridiagonal symmetric system Ax = b
% A has diagonal d and subdiagonal beta
% Cost: O(m1) — forward sweep + back substitution
%
% INPUT:
%   d    - modified diagonal (m1 x 1)
%   beta - subdiagonal (m1-1 x 1)
%   b    - right-hand side (m1 x 1)
%
% OUTPUT:
%   x    - solution (m1 x 1)

    m1 = length(d);

    %  Forward sweep (LU decomposition) 
    % Modify diagonal in-place
    d_mod = d;       % copy to avoid modifying the original
    b_mod = b;       % copy of right-hand side

    for j = 2:m1
        % Elimination factor
        factor = beta(j-1) / d_mod(j-1);
        % Update diagonal
        d_mod(j) = d_mod(j) - factor * beta(j-1);
        % Update right-hand side
        b_mod(j) = b_mod(j) - factor * b_mod(j-1);
    end

    %  Back substitution 
    x = zeros(m1, 1);
    x(m1) = b_mod(m1) / d_mod(m1);

    for j = m1-1:-1:1
        x(j) = (b_mod(j) - beta(j) * x(j+1)) / d_mod(j);
    end

end