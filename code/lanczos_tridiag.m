function [Phi, alpha, beta, m1] = lanczos_tridiag(R, rho, l_k)
% LANCZOS_TRIDIAG  Lanczos tridiagonalization of R

% INPUT:
%   R    - correlation matrix (m x m), symmetric positive semidefinite
%   rho  - stopping threshold (default 0.995)
%   l_k  - number of principal components needed
%
% OUTPUT:
%   Phi   - Lanczos basis vectors (m x m1), orthonormal columns
%   alpha - diagonal of tridiagonal Gamma (m1 x 1)
%   beta  - subdiagonal of Gamma (m1-1 x 1)
%   m1    - number of Lanczos iterations performed


    m = size(R, 1);           % number of variables
    tr_R = trace(R);          % = m for normalized data (fixed, no eig needed)
    
    % Maximum iterations: never exceed m
    % Minimum iterations: at least l_k (need l_k eigenvalues)
    max_iter = m;
    
    % Pre-allocate (maximum possible size)
    Phi   = zeros(m, max_iter);  % Lanczos vectors
    alpha = zeros(max_iter, 1);  % diagonal
    beta  = zeros(max_iter, 1);  % subdiagonal (beta(j) = beta_j)
    
    % Step 1: Initial vector 
    % Random unit vector
    q = randn(m, 1);
    q = q / norm(q);
    Phi(:, 1) = q;
    
    % Step 2: Lanczos iteration
    tr_Gamma = 0;  % accumulate trace of Gamma
    m1 = 0;        % iteration counter
    
    for j = 1:max_iter
        
        m1 = j;
        
        % Current Lanczos vector
        q_j = Phi(:, j);
        
        % Compute R * q_j
        r = R * q_j;
        
        % Compute alpha_j = q_j' * R * q_j  (Rayleigh quotient)
        alpha(j) = q_j' * r;
        
        % Subtract alpha_j * q_j component
        r = r - alpha(j) * q_j;
        
        % Subtract beta_{j-1} * q_{j-1} component (if j > 1)
        if j > 1
            r = r - beta(j-1) * Phi(:, j-1);
        end
        
        % Re-orthogonalization (full) 
        % Numerically necessary: loss of orthogonality in floating point
        % Project out all previous Lanczos vectors
        for i = 1:j
            r = r - (Phi(:,i)' * r) * Phi(:,i);
        end
        
        % Compute beta_j = ||r||
        beta(j) = norm(r);
        
        % Check stopping criterion
        % Update trace of Gamma
        tr_Gamma = tr_Gamma + alpha(j);
        
        % Stop if:
        % (1) captured enough variance
        % (2) OR beta_j is numerically zero (R-invariant subspace reached)
        % (3) OR we have at least l_k iterations
        
        if (tr_Gamma / tr_R >= rho) && (j >= l_k)
            break
        end
        
        if beta(j) < 1e-10
            % R-invariant subspace: Krylov space exhausted
            fprintf('Lanczos: R-invariant subspace at j=%d\n', j);
            break
        end
        
        % Next Lanczos vector
        if j < max_iter
            Phi(:, j+1) = r / beta(j);
        end
        
    end
    
    % Trim outputs to actual size m1
    Phi  = Phi(:, 1:m1);
    alpha = alpha(1:m1);
    beta  = beta(1:m1-1);  % beta has m1-1 elements
    
end