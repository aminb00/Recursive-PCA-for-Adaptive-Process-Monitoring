function c = inverse_iteration(alpha, beta, lambda_i, delta)
% INVERSE_ITERATION  Trova l'autovettore di Gamma per lambda_i
% usando inverse iteration con Thomas algorithm
% Following Li et al. (2000), Section 4
%
% INPUT:
%   alpha    - diagonale di Gamma (m1 x 1)
%   beta     - sottodiagonale di Gamma (m1-1 x 1)
%   lambda_i - autovalore noto (da bisection_eig)
%   delta    - perturbazione per evitare singolarità (default 1e-10)
%
% OUTPUT:
%   c - autovettore normalizzato (m1 x 1)

    if nargin < 4
        delta = 1e-10;
    end

    m1 = length(alpha);

    % --- Matrice shiftata: A = Gamma - lambda_i*I + delta*I ---
    % Diagonale shiftata
    d = alpha - lambda_i + delta;  % (m1 x 1)
    % Sottodiagonale rimane beta

    % --- Vettore iniziale casuale normalizzato ---
    c = randn(m1, 1);
    c = c / norm(c);

    % --- Inverse iteration: 3 iterazioni sono sempre sufficienti ---
    % se lambda_i è noto con precisione 1e-10
    max_iter = 3;

    for iter = 1:max_iter

        % Risolvi (Gamma - lambda_i*I + delta*I) * c_new = c
        % cioè: A * c_new = c
        % con A tridiagonale simmetrica → Thomas algorithm

        c_new = thomas_tridiag(d, beta, c);

        % Normalizza
        c_new_norm = norm(c_new);

        if c_new_norm < 1e-14
            % Numericamente zero — problema di convergenza
            warning('inverse_iteration: norma quasi zero a iter %d', iter);
            break
        end

        c = c_new / c_new_norm;

    end

end


function x = thomas_tridiag(d, beta, b)
% THOMAS_TRIDIAG  Risolve sistema tridiagonale simmetrico Ax = b
% A ha diagonale d e sottodiagonale beta
% Costo: O(m1) — forward sweep + back substitution
%
% INPUT:
%   d    - diagonale modificata (m1 x 1)
%   beta - sottodiagonale (m1-1 x 1)
%   b    - termine noto (m1 x 1)
%
% OUTPUT:
%   x    - soluzione (m1 x 1)

    m1 = length(d);

    % --- Forward sweep (LU decomposition) ---
    % Modifica diagonale in-place
    d_mod = d;       % copia per non modificare l'originale
    b_mod = b;       % copia del termine noto

    for j = 2:m1
        % Fattore di eliminazione
        factor = beta(j-1) / d_mod(j-1);
        % Aggiorna diagonale
        d_mod(j) = d_mod(j) - factor * beta(j-1);
        % Aggiorna termine noto
        b_mod(j) = b_mod(j) - factor * b_mod(j-1);
    end

    % --- Back substitution ---
    x = zeros(m1, 1);
    x(m1) = b_mod(m1) / d_mod(m1);

    for j = m1-1:-1:1
        x(j) = (b_mod(j) - beta(j) * x(j+1)) / d_mod(j);
    end

end